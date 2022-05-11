#ifndef RK45_HPP
#define RK45_HPP

#include <iostream>
#include <vector>   // for std::vector
#include <utility>  // for std::pair

#include "vector.hpp"

// A class to hold the ODE solver with the option to define events during integration.
namespace integrator
{
    typedef vector (*ODE_system)(const double, const vector& , const void*);
    typedef std::pair<double, vector> step;
    typedef bool (*event_condition)(const double, const vector, const void* params);

    // holds an "Event" which can be checked during integration (e.g. a stopping condition)
    struct Event {
        event_condition condition;
        bool stopping_condition;
        std::vector<step> steps;
        bool active;
        std::string name;
        Event(event_condition condition, bool stopping_condition=false, std::string name="") : condition(condition), stopping_condition(stopping_condition),
                                                                        active(false), name(name) {}
        void reset() { steps.clear(); active=false; }
    };

    // holds parameters/options relevant to the integrator
    struct IntegrationOptions {
        double target_error;
        int max_step;
        double min_stepsize;
        double max_stepsize;
        bool save_intermediate;
        IntegrationOptions(const int max_step=1000000, const double target_error=1e-10, const double min_stepsize=1e-18, const double max_stepsize=1., const bool save_intermediate=false)
                            : max_step(max_step), target_error(target_error), min_stepsize(min_stepsize), max_stepsize(max_stepsize), save_intermediate(save_intermediate) {}
    };

    // define return codes fpr the integrator
    enum  return_reason  {endpoint_reached=1, stepsize_underflow, iteration_number_exceeded,  event_stopping_condition};

    // Runge-Kutta Fehlberg stepper
    bool RKF45_step(ODE_system dy_dt, double &r, double &dr, vector& y, const void* params, const IntegrationOptions& options);
    // Full Runge-Kutta Fehlberg IVP integrator
    int RKF45(ODE_system dy_dt, const double r0, const vector y0, const double r_end, const void* params,
                            std::vector<step>& results, std::vector<Event>& events, const IntegrationOptions& options);

}

#endif