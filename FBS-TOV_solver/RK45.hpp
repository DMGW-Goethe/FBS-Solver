#ifndef RK45_HPP
#define RK45_HPP

#include <iostream>
#include <vector>
#include <utility>

#include "vec5d.hpp"
#include "eostable.hpp"


namespace RK45
{
    typedef vec5d vec;
    typedef vec (*ODE_system)(const double, const vec& , const void*);
    typedef std::pair<double, vec> step;
    typedef bool (*event_condition)(const double, const vec, const void* params);

    void RKF45_step(ODE_system dy_dt, double &r, double &dr, vec& y, const void* params, const double target_error=1e-9, const double max_stepsize=1e-3, const double min_stepsize=1e-5);
    int RKF45(ODE_system dy_dt, const double r0, const vec y0, const double r_end, const void* params, const int max_step,
                            std::vector<step>& results, std::vector<step>& events, const bool save_intermediate=false, event_condition event=NULL, event_condition stopping_condition=NULL);

    void RK45_Fehlberg(ODE_system dy_dt, double &r, double &dr , vec5d& V, const void *params);

}





#endif