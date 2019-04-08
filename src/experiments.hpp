//
// Created by mitja on 22.3.2019.
//

#ifndef VAJA_II_2_EXPERIMENTS_HPP
#define VAJA_II_2_EXPERIMENTS_HPP

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <algorithm>
#include "SimplecticS4.h"

static long seed = std::chrono::system_clock::now().time_since_epoch().count();
static std::default_random_engine generator (seed);
static std::normal_distribution<double> distribution (0.0,1.0);

static const size_t N = 10;
static const size_t DIM = 2*N + 2;
static const double binDelimiter = -1234567891.0;

void writeBinary(std::vector<double>& data, const std::string& file)
{
    std::ofstream output(file, std::ios::binary);
    for (double& x : data)
    {
        output.write(reinterpret_cast<char*>(&x), sizeof(x));
    }
    output.close();
}

int systemFunc(double t, const double y[], double f[], void * params)
{
    (void) (t);
    double relax  = *(double*)params;
    double tl = *((double*)params + 1);
    double tr = *((double*)params + 2);
    double lambda = *((double*)params + 3);

    for (size_t j = 0; j < N; ++j)
    {
        f[j] = y[j + N + 1];
    }

    f[N] = relax * (y[N + 1] * y[N + 1] - tl);
    f[N + 1] = -2 * y[0] - 4 * lambda * pow(y[0], 3) + y[1] - y[N] * y[N + 1];

    for (size_t j = 1; j < N - 1; ++j)
    {
        f[j + N + 1] = -3 * y[j] - 4 * lambda * pow(y[j], 3) + y[j - 1] + y[j + 1];
    }

    f[2*N] = -2 * y[N - 1] - 4 * lambda * pow(y[N - 1], 3) + y[N - 2] - y[2*N + 1] * y[2*N];
    f[2*N + 1] = relax * (y[2*N] * y[2*N] - tr);

    return GSL_SUCCESS;
}

void stateInit(double y[DIM])
{
    for (size_t j = 0; j < DIM; ++j)
    {
        if (j > N && j <= 2*N) y[j] = distribution(generator);
        else y[j] = 0;
    }
}

void stateInit(std::vector<double>& x)
{
    for (size_t j = N; j < x.size(); ++j)
    {
        x[j] = distribution(generator);
    }
}

void addTemperatures (const double y[DIM], std::vector<double>& target)
{
    for (size_t j = 0; j < N; ++j)
    {
        target[j] += y[j + N + 1] * y[j + N + 1] / 2;
    }
}

void makeTProfile(double step, size_t steps, double *params, double *y, std::vector<double>& target)
{
    gsl_odeiv2_system sys = {systemFunc, nullptr, DIM, params};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 0.0, 1e-3);

    double t1 = 0;

    std::vector<double> accumulate(N, 0);

    for (int t = 0; t < steps; ++t)
    {
        std::cout << "Step #" << t << std::endl;
        int status = gsl_odeiv2_driver_apply (d, &t1, t1 + step, y);
        if (status != GSL_SUCCESS)
        {
            std::cerr << "error, return value=" << status << std::endl;
            break;
        }
        target.push_back(t1);
        addTemperatures(y, accumulate);
        std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                [t1, step](double x){return step*x/t1;});
        target.push_back(binDelimiter);
    }
}

void addFluxes (const double y[DIM], std::vector<double>& target)
{
    for (size_t j = 1; j < N-1; ++j)
    {
        target[j - 1] += (y[j - 1] - y[j + 1]) * y[j + N + 1] / 2;
    }
}

void makeFlux(double step, size_t steps, double *params, double *y, std::vector<double>& target)
{
    gsl_odeiv2_system sys = {systemFunc, nullptr, DIM, params};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 0.0, 1e-3);

    double t1 = 0;

    std::vector<double> accumulate(N-2, 0);

    for (int t = 0; t < steps; ++t)
    {
        std::cout << "Step #" << t << std::endl;
        int status = gsl_odeiv2_driver_apply (d, &t1, t1 + step, y);
        if (status != GSL_SUCCESS)
        {
            std::cerr << "error, return value=" << status << std::endl;
            break;
        }
        target.push_back(t1);
        addFluxes(y, accumulate);
        std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                [t1, step](double x){return step*x/t1;});
        target.push_back(binDelimiter);
    }
}

void maxwelTProfile(double tau, size_t reads, size_t samples, size_t steps, double lambda,
                    std::normal_distribution<double>& tl, std::normal_distribution<double>& tr,
                    std::vector<double>& target)
{
    SimplecticS4 integrator(tau / steps, lambda);
    std::vector<double> x(2*N, 0);
    stateInit(x);

    double time = 0;
    double readStep = tau * samples;
    std::vector<double> accumulate(N, 0);

    for (size_t read = 0; read < reads; ++read)
    {
        std::cout << "Read #" << read << std::endl;

        for (size_t sample = 0; sample < samples; ++sample)
        {
            integrator.propagate(x, steps);

            x[N] = tl(generator);
            x[2*N - 1] = tr(generator);
        }

        time += readStep;
        target.push_back(time);

        std::transform(accumulate.begin(), accumulate.end(), x.begin() + N, accumulate.begin(),
                       [](double acc, double x){return acc + x*x/2;});
        std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                       [time, readStep](double x){return readStep*x/time;});

        target.push_back(binDelimiter);
    }
}


#endif //VAJA_II_2_EXPERIMENTS_HPP
