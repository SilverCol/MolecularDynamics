//
// Created by mitja on 22.3.2019.
//

#ifndef VAJA_II_2_EXPERIMENTS_HPP
#define VAJA_II_2_EXPERIMENTS_HPP

#include <iostream>
#include <chrono>
#include <fstream>
#include <cmath>

static const size_t N = 50;

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
}

#endif //VAJA_II_2_EXPERIMENTS_HPP
