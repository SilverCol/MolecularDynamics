//
// Created by mitja on 22.3.2019.
//

#ifndef VAJA_II_2_EXPERIMENTS_HPP
#define VAJA_II_2_EXPERIMENTS_HPP

#include <iostream>
#include <chrono>
#include <fstream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

static const size_t N = 50;
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
    for (size_t j = 0; j < DIM; ++j) y[j] = 0;
}

void addTemperatures(const double y[DIM], std::vector<double>& target)
{
    for (size_t j = N + 1; j < 2*N + 1; ++j)
    {
        target.push_back(y[j] * y[j] / 2);
    }
}

void makeTProfile(double step, size_t steps, double params[4], double y[DIM], std::vector<double>& target)
{
    gsl_odeiv2_system sys = {systemFunc, nullptr, DIM, params};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 0.0, 1e-6);

    double t1 = 0;
    double t2 = step;

    target.push_back(t1);
    addTemperatures(y, target);
    target.push_back(binDelimiter);

    for (int t = 0; t < steps; ++t)
    {
        int status = gsl_odeiv2_driver_apply (d, &t1, t2, y);
        if (status != GSL_SUCCESS)
        {
            std::cerr << "error, return value=" << status << std::endl;
            break;
        }
        t2 += step;
        target.push_back(t1);
        addTemperatures(y, target);
        target.push_back(binDelimiter);
    }
}

#endif //VAJA_II_2_EXPERIMENTS_HPP
