//
// Created by mitja on 7.4.2019.
//

#include <cmath>
#include "SimplecticS4.h"

SimplecticS4::SimplecticS4(double tau, double lambda):
m_tau(tau),
m_lambda(lambda),
m_x0(-1.7024143839193153 * tau),
m_x1(1.3512071919596578 * tau),
m_x1h(m_x1/2),
m_x01((m_x0 + m_x1)/2)
{}

void SimplecticS4::propagate(std::vector<double>& x, size_t steps)
{
    for (size_t n = 0; n < steps; ++n)
    {
        step(x);
    }
}

void SimplecticS4::step(std::vector<double>& x)
{
    size_t N = x.size() / 2;

    kinetic(x, N, m_x1h);
    potential(x, N, m_x1);
    kinetic(x, N, m_x01);
    potential(x, N, m_x0);
    kinetic(x, N, m_x01);
    potential(x, N, m_x1);
    kinetic(x, N, m_x1h);
}

void SimplecticS4::kinetic(std::vector<double>& x, size_t N, double factor)
{
    for (size_t j = 0; j < N; ++j) x[j] += factor * x[j + N];
}

void SimplecticS4::potential(std::vector<double>& x, size_t N, double factor)
{
    x[N] -= factor * (2*x[0] + 4 * m_lambda * pow(x[0], 3) - x[1]);
    for (size_t j = 1; j < N - 1; ++j) x[N + j] -= factor * (2*x[j] + 4 * m_lambda * pow(x[j], 3) - x[j - 1] - x[j + 1]);
    x[2*N - 1] -= factor * (2*x[N - 1] + 4 * m_lambda * pow(x[N - 1], 3) - x[N - 2]);
}
