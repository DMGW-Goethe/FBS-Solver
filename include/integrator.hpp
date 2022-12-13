#ifndef RK45_HPP
#define RK45_HPP

#include <iostream>
#include <vector>   // for std::vector
#include <utility>  // for std::pair
#include <stdexcept> // for std::runtime_error

#include "vector.hpp"

// A class to hold the ODE solver with the option to define events during integration.
namespace integrator
{
    /* Define the ODE_system value, which is e.g. the dy_dr in the NSmodel class */
    typedef vector (*ODE_system)(const double, const vector& , const void*);

    /* the integrator works in steps of r and the integrated variables */
    typedef std::pair<double, vector> step;

    /* function pointer to event conditions, tests whether a condition is fulfilled */
    typedef bool (*event_condition)(const double r, const double dr, const vector& y, const vector& dy, const void*params);

    /* Event
     * describes an integration "event" which can be checked during integration (e.g. a stopping condition)
     * has the following parameters
     *  condition           : pointer to a function checking whether the condition is fulfilled
     *  stopping_condition  : a bool whether this is a stopping condition, i.e. the integrator will stop once this is fulfilled
     *                          in that case event_stopping_condition is return_reason, then the causing event can be found out with the active field
     *  steps               : A list of steps where the event BECAME active (NOT all points where event is active)
     *  name                : An optional name parameter
     *  active              : Is a field used by the integrator to save whether it is active or not */
    struct Event {
        event_condition condition;
        bool stopping_condition;
        std::vector<step> steps;
        bool active;
        std::string name;
        double target_accuracy;
        Event(event_condition condition, bool stopping_condition=false, std::string name="", double target_accuracy=0.) : condition(condition), stopping_condition(stopping_condition),
                                                                        active(false), name(name), target_accuracy(target_accuracy) {}
        /* This resets the event before an integration - whether this is called depends on the clean_events parameter in the IntegrationOptions */
        void reset() { steps.clear(); active=false; }
    };

    /* Simple output for the step type */
    std::ostream& operator <<(std::ostream& o, const step& s);

    /* IntegrationOption
     * contains options for the integration
     *  max_step      : the maximal number of steps the integrator should take - if exceeded iteration_number_exceeded is the return_reason
     *  min_stepsize  : the minimal stepsize for an integration step - if smaller steps are needed an stepsize_underflow is the return_reason
     *  max_stepsize  : the maximal stepsize for an integration step
     *  save_intermediate : Whether the integrator should save all steps in the results vector or just the initial and last steps
     *  verbose       : How verbose the integrator should be
     *  clean_events  : Whether to call event->reset before integration - Default is true, put to false when continuing an integration for example
     */
    struct IntegrationOptions {
        int max_step;
        double target_error;
        double min_stepsize;
        double max_stepsize;
        bool save_intermediate;
        int verbose;
        bool clean_events;
        IntegrationOptions(const int max_step=1000000, const double target_error=1e-12, const double min_stepsize=1e-16, const double max_stepsize=1e-2, const bool save_intermediate=false, const int verbose=0, bool clean_events=true)
                            : max_step(max_step), target_error(target_error), min_stepsize(min_stepsize), max_stepsize(max_stepsize), save_intermediate(save_intermediate), verbose(verbose), clean_events(clean_events) {}
    };

    /* define return codes for the integrator */
    enum  return_reason  {endpoint_reached=1, stepsize_underflow, iteration_number_exceeded,  event_stopping_condition};

    /* Runge-Kutta Fehlberg stepper that does one step at a time */
    bool RKF45_step(ODE_system dy_dr, double &r, double &dr, vector& y, const void* params, const IntegrationOptions& options);

    /* Checks if events require smaller stepsize and only accepts steps if they do*/
    int RKF45_step_event_tester(ODE_system dy_dr, step& current_step, double& dr, const void* params,
                                            const std::vector<Event>& events, const IntegrationOptions& options);

    /* Full Runge-Kutta Fehlberg IVP integrator, steps are output in results vector*/
    int RKF45(ODE_system dy_dr, const double r0, const vector y0, const double r_end, const void* params,
                            std::vector<step>& results, std::vector<Event>& events, const IntegrationOptions& options);

    /* A simple cumulative trapezoid integration algorithm */
    void cumtrapz(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& res);

}

#endif
