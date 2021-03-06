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
static std::normal_distribution<double> distribution (0.0,1.225);

static const size_t N = 80;
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

void stateInit(std::vector<double>& x)
{
    for (size_t j = x.size()/2; j < x.size()-1; ++j)
    {
        x[j] = distribution(generator);
    }
}

void maxwellInit(std::vector<double>& x)
{
    for (size_t j = x.size()/2; j < x.size(); ++j)
    {
        x[j] = distribution(generator);
    }
}

void makeProfile(double step, size_t steps, double *params, std::vector<double>& target, std::vector<double>& x,
                  double cutoff)
{
    gsl_odeiv2_system sys = {systemFunc, nullptr, x.size(), params};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 0.0, 1e-3);

    double t = 0;
    std::cout << "Cutting the front end: " << cutoff << std::endl;
    int status = gsl_odeiv2_driver_apply (d, &t, cutoff, x.data());
    if (status != GSL_SUCCESS)
    {
        std::cerr << "error, return value=" << status << std::endl;
    }
    t = 0;

    std::vector<double> accumulate(N, 0);

    for (size_t i = 0; i < steps; ++i)
    {
        std::cout << "Step #" << i << std::endl;
        int status = gsl_odeiv2_driver_apply (d, &t, t + step, x.data());
        if (status != GSL_SUCCESS)
        {
            std::cerr << "error, return value=" << status << std::endl;
            break;
        }
        std::transform(accumulate.begin(), accumulate.end(), x.begin() + N + 1, accumulate.begin(),
                       [](double acc, double x){return acc + x*x;});
    }

    std::cout << "Elapsed time: " << t << std::endl;
    std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                   [t, step](double x){return step*x/t;});
}

void makeTProfile(double step, size_t steps, double *params, std::vector<double>& target, std::vector<double>& x,
                  double cutoff)
{
    gsl_odeiv2_system sys = {systemFunc, nullptr, x.size(), params};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 0.0, 1e-3);

    double t = 0;
    std::cout << "Cutting the front end: " << cutoff << std::endl;
    int status = gsl_odeiv2_driver_apply (d, &t, cutoff, x.data());
    if (status != GSL_SUCCESS)
    {
        std::cerr << "error, return value=" << status << std::endl;
    }
    t = 0;

    std::vector<double> accumulate(N, 0);

    for (size_t i = 0; i < steps; ++i)
    {
        std::cout << "Step #" << i << std::endl;
        int status = gsl_odeiv2_driver_apply (d, &t, t + step, x.data());
        if (status != GSL_SUCCESS)
        {
            std::cerr << "error, return value=" << status << std::endl;
            break;
        }

        target.push_back(t);
        std::transform(accumulate.begin(), accumulate.end(), x.begin() + N + 1, accumulate.begin(),
                       [](double acc, double x){return acc + x*x;});
        std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                [t, step](double x){return step*x/t;});
        target.push_back(binDelimiter);
    }
}

void addFluxes (const std::vector<double>& y, std::vector<double>& target)
{
    for (size_t j = 1; j < N-1; ++j)
    {
        target[j - 1] += (y[j - 1] - y[j + 1]) * y[j + N + 1] / 2;
    }
}

void maxwellFluxes (const std::vector<double>& y, std::vector<double>& target)
{
    for (size_t j = 1; j < N-1; ++j)
    {
        target[j - 1] += (y[j - 1] - y[j + 1]) * y[j + N] / 2;
    }
}

void flux(double step, size_t steps, double *params, std::vector<double>& target, std::vector<double>& x)
{
    gsl_odeiv2_system sys = {systemFunc, nullptr, x.size(), params};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 0.0, 1e-2);

    double t = 0;

    std::vector<double> accumulate(N-2, 0);

    for (size_t i = 0; i < steps; ++i)
    {
        std::cout << "Step #" << i << std::endl;
        int status = gsl_odeiv2_driver_apply (d, &t, t + step, x.data());
        if (status != GSL_SUCCESS)
        {
            std::cerr << "error, return value=" << status << std::endl;
            break;
        }
        addFluxes(x, accumulate);
    }
    std::cout << "Elapsed time: " << t << std::endl;
    std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                   [t, step](double x){return step*x/t;});
}

void makeFlux(double step, size_t steps, double *params, std::vector<double>& target, std::vector<double>& x)
{
    gsl_odeiv2_system sys = {systemFunc, nullptr, x.size(), params};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 0.0, 1e-3);

    double t = 0;

    std::vector<double> accumulate(N-2, 0);

    for (size_t i = 0; i < steps; ++i)
    {
        std::cout << "Step #" << i << std::endl;
        int status = gsl_odeiv2_driver_apply (d, &t, t + step, x.data());
        if (status != GSL_SUCCESS)
        {
            std::cerr << "error, return value=" << status << std::endl;
            break;
        }
        target.push_back(t);
        addFluxes(x, accumulate);
        std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                [t, step](double x){return step*x/t;});
        target.push_back(binDelimiter);
    }
}

void maxwelProfile(double tau, size_t reads, size_t samples, size_t steps, double lambda,
                    std::normal_distribution<double>& tl, std::normal_distribution<double>& tr,
                    std::vector<double>& target, std::vector<double>& x, double cutoff)
{
    SimplecticS4 integrator(tau / steps, lambda);

    double time = 0;
    std::cout << "Cutting the front end: " << cutoff << std::endl;
    integrator.propagate(x, steps*(size_t)(cutoff/tau));

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
        std::transform(accumulate.begin(), accumulate.end(), x.begin() + N, accumulate.begin(),
                       [](double acc, double x){return acc + x*x;});
    }
    std::cout << "Elapsed time: " << time << std::endl;
    std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                   [time, readStep](double x){return readStep*x/time;});
}

void maxwelTProfile(double tau, size_t reads, size_t samples, size_t steps, double lambda,
                    std::normal_distribution<double>& tl, std::normal_distribution<double>& tr,
                    std::vector<double>& target, std::vector<double>& x, double cutoff)
{
    SimplecticS4 integrator(tau / steps, lambda);

    double time = 0;
    std::cout << "Cutting the front end: " << cutoff << std::endl;
    integrator.propagate(x, steps*(size_t)(cutoff/tau));

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
                       [](double acc, double x){return acc + x*x;});
        std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                       [time, readStep](double x){return readStep*x/time;});

        target.push_back(binDelimiter);
    }
}

void fluxMax(double tau, size_t reads, size_t samples, size_t steps, double lambda,
                    std::normal_distribution<double>& tl, std::normal_distribution<double>& tr,
                    std::vector<double>& target, std::vector<double>& x)
{
    SimplecticS4 integrator(tau / steps, lambda);

    double time = 0;
    double readStep = tau * samples;
    std::vector<double> accumulate(N - 2, 0);

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
        maxwellFluxes(x, accumulate);
    }
    std::cout << "Elapsed time: " << time << std::endl;
    std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                   [time, readStep](double x){return readStep*x/time;});
}

void maxwelFlux(double tau, size_t reads, size_t samples, size_t steps, double lambda,
                    std::normal_distribution<double>& tl, std::normal_distribution<double>& tr,
                    std::vector<double>& target, std::vector<double>& x)
{
    SimplecticS4 integrator(tau / steps, lambda);

    double time = 0;
    double readStep = tau * samples;
    std::vector<double> accumulate(N - 2, 0);

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

        maxwellFluxes(x, accumulate);
        std::transform(accumulate.begin(), accumulate.end(), std::back_inserter(target),
                       [time, readStep](double x){return readStep*x/time;});

        target.push_back(binDelimiter);
    }
}

void contactFlux(double itau, double dtau, size_t nTau, size_t samples, size_t steps, double lambda,
                    std::normal_distribution<double>& tl, std::normal_distribution<double>& tr,
                    std::vector<double>& target, std::vector<double>& x)
{
    double tau = itau;
    for (size_t i = 0; i < nTau; ++i)
    {
        SimplecticS4 integrator(tau / steps, lambda);

        std::cout << "Tau = " << tau << std::endl;

        double j1 = 0;
        double j2 = 0;
        for (size_t sample = 0; sample < samples; ++sample)
        {
            integrator.propagate(x, steps);

            j1 -= x[N]*x[N];
            x[N] = tl(generator);
            j1 += x[N]*x[N];

            j2 -= x[2*N - 1]*x[2*N - 1];
            x[2 * N - 1] = tr(generator);
            j2 += x[2*N - 1]*x[2*N - 1];
        }

        target.push_back(tau);
        target.push_back(j1/(2*samples*tau));
        target.push_back(j2/(2*samples*tau));
        tau += dtau;
    }
}

#endif //VAJA_II_2_EXPERIMENTS_HPP
