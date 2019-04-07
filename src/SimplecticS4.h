//
// Created by mitja on 7.4.2019.
//

#ifndef VAJA_III_1_SIMPLECTICS4_H
#define VAJA_III_1_SIMPLECTICS4_H

#include <vector>

class SimplecticS4
{
    SimplecticS4(double tau, double lambda);
    void propagate(std::vector<double>& x, size_t steps);
private:
    void step(std::vector<double>& x);
    void kinetic(std::vector<double>& x, size_t N, double factor);
    void potential(std::vector<double>& x, size_t N, double factor);
    double m_tau;
    double m_lambda;
    double m_x0;
    double m_x1;
    double m_x1h;
    double m_x01;
};


#endif //VAJA_III_1_SIMPLECTICS4_H
