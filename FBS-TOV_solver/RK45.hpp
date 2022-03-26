#ifndef RK45_HPP
#define RK45_HPP

#include <iostream>
#include <vector>
#include <utility>

#include "vector.hpp"


namespace integrator
{
    typedef vector (*ODE_system)(const double, const vector& , const void*);
    typedef std::pair<double, vector> step;
    typedef bool (*event_condition)(const double, const vector, const void* params);

    struct Event {
        event_condition condition;
        bool stopping_condition;
        std::vector<step> steps;
        bool active;
        std::string name;
        Event(event_condition condition, bool stopping_condition=false, std::string name="") : condition(condition), stopping_condition(stopping_condition),
                                                                        active(false), name(name) {}
        void reset() { steps.empty(); active=false; }
    };

    void RKF45_step(ODE_system dy_dt, double &r, double &dr, vector& y, const void* params, const double target_error=1e-9, const double max_stepsize=1e-3, const double min_stepsize=1e-5);
    int RKF45(ODE_system dy_dt, const double r0, const vector y0, const double r_end, const void* params, const int max_step,
                            std::vector<step>& results, std::vector<Event>& events, const bool save_intermediate=false);

}





#endif
