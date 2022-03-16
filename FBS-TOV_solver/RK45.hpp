#ifndef RK45_HPP
#define RK45_HPP

#include <iostream>
#include <vector>

#include "vec5d.hpp"
#include "eostable.hpp"


namespace RK45
{
    typedef vec5d vec;
    typedef vec (*ODE_system)(const double, const vec& , const void*);
//typedef bool (*stopping_condition)(const double, const vec, const void* params);

//void RKF45_step(ode_system ode, double &r, double &dr, vec& y, const void* params);
//void RKF45(ode_system ode, const double r0, const vec y0, const double r_end, const void* params, const int max_step);

//void RK45_Fehlberg(double &r, double &dr , vec5d& V, eostable& myEOS, const double mu, const double lambda, const double omega);
    void RK45_Fehlberg(ODE_system dy_dt, double &r, double &dr , vec5d& V, const void *params);

}





#endif
