#include "integrator.hpp"

namespace ublas = boost::numeric::ublas;

/* This function takes as input the ODE_system, the position r, the current stepsize dr, and the parameters y
 * It perfoms one RKF45 step and saves the results in the updated r, stepsize dr and y
 * If the stepsize is large enough it is reduced for the next step
 * */
bool integrator::RKF45_step(ODE_system dy_dr, NUMERIC &r, NUMERIC &dr, vector& y, const void* params, const IntegrationOptions& options)
{
    // intermediate steps:
    int n = y.size();
    vector k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), dV_04(n), dV_05(n);

    // solve until we find an acceptable stepsize/error:
    while (true) {

        // runge kutta Fehlberg steps for all ODEs:
        k1 = dr * dy_dr(r, y, params);
        k2 = dr * dy_dr(r + 1._num / 4._num * dr, y + 1._num / 4._num * k1, params);
        k3 = dr * dy_dr(r + 3._num / 8._num * dr, y + 3._num / 32._num * k1 + 9._num / 32._num * k2, params);
        k4 = dr * dy_dr(r + 12._num / 13._num * dr, y + 1932._num / 2197._num * k1 - 7200._num / 2197._num * k2 + 7296._num / 2197._num * k3, params);
        k5 = dr * dy_dr(r + dr, y + 439._num / 216._num * k1 - 8._num * k2 + 3680._num / 513._num * k3 - 845._num / 4104._num * k4, params);
        k6 = dr * dy_dr(r + 1._num / 2._num * dr, y - 8._num / 27._num * k1 + 2._num * k2 - 3544._num / 2565._num * k3 + 1859._num / 4104._num * k4 - 11._num / 40._num * k5, params);

        // 4th and 5th order accurate steps:
		dV_04 = 25._num / 216._num * k1 + 1408._num / 2565._num * k3 + 2197._num / 4104._num * k4 - 1._num / 5._num * k5; // + O(x^5)
        dV_05 = 16._num / 135._num * k1 + 6656._num / 12825._num * k3 + 28561._num / 56430._num * k4 - 9._num / 50._num * k5 + 2._num / 55._num * k6; // + O(x^6)

        if( vector::is_nan(dV_05) || vector::is_nan(dV_04)) {
            dr *= 0.5_num;
            if (dr < options.min_stepsize) {
                if(options.verbose > 0)
                    std::cout << "Nan found" << dV_05 << dV_04 << k1 << k2 << k3 << k4 << k5 << k6 << "\n  for r=" << r << ", dr=" << dr <<  std::endl;
                throw std::runtime_error("NaN detected");
            }
            if(options.verbose > 3)
                std::cout << "step not precise enough and stepsize decreased creased: dr = " << dr << std::endl;
            continue;
        }

		if (options.force_max_stepsize) {	// note: this will ignore target truncation error
			r += dr;
            y += dV_05;
			dr = options.max_stepsize;
			return true;
		}

        // approximating the truncation error:
        NUMERIC truncation_error = ublas::norm_inf(dV_05 - dV_04) * dr; // inf-Norm
		//NUMERIC truncation_error = ublas::norm_2(dV_05 - dV_04) * dr; // 2-Norm

		if (truncation_error > options.target_error) {
			dr *= 0.5_num;     // truncation error is too large. We repeat the iteration with smaller stepsize

            if (dr < options.min_stepsize) {
                // error is not acceptable but the stepsize cannot get any smaller:
                // write new values for r and the evolved quantities
                r += dr;
                y += dV_05;

                dr = options.min_stepsize; // ensure that the stepsize never gets too small
                if(options.verbose > 0) {
                    std::cout << "Error in RKF45_step(): Minimal stepsize underflow at r = " << r << ", dr = "<< dr << std::endl;  // error message for debugging
                    std::cout << "Truncation error = "<< truncation_error << std::endl;
                    std::cout << y << std::endl;
                }
                return false;
            }

            if(options.verbose > 3)
                std::cout << "step not precise enough and stepsize decreased creased: dr = " << dr << std::endl;
			continue;
		}
		else {
			// error is small enough and therefore we can use this step. To save computation time, we then increase the stepsize
            r += dr;
            y += dV_05;

            dr *= 2.0_num;
			if (dr > options.max_stepsize)
                dr = options.max_stepsize;   // enforce maximal stepsize

            if(options.verbose > 3){
                std::cout << "step accepted and stepsize increased: dr = " << dr << std::endl;
            }
			return true;
		}

	    return true;  // stop this iteration step if the required accuracy and step size was acheived
	}
}

int integrator::RKF45_step_event_tester(ODE_system dy_dr, step& current_step, NUMERIC& step_size, const void* params,
                                            const std::vector<Event>& events, const IntegrationOptions& options) {

    NUMERIC r;
    NUMERIC dr;
    vector y, dy;
    bool step_success = false;

    while(!step_success ) {
        r = current_step.first;
        y = current_step.second;
        step_success = RKF45_step(dy_dr, r, step_size, y, params, options);
        if (!step_success) // if this step fails already, return
            break;

        dr = r - current_step.first;
        dy = dy_dr(r, y, params);
        // iterate through all defined events and check if they would turn active
        for(auto it = events.begin(); it != events.end(); ++it) {
            if ( it->active || it->target_accuracy <= 0._num)
                continue;

            if(it->condition(r, dr, y, dy, params)) {
                if ( dr > it->target_accuracy*1.001_num ) { // do these require higher accuracy?
                    step_success = false;
                    if (options.verbose > 1)
                        std::cout << "event " << it->name << " required higher accuracy at r = " << r << " where stepsize = " << dr << " but req " << it->target_accuracy << std::endl;
                }
            }
        }
        if (step_success) // none would be active or target accuracy is achieved
            break;

        if (dr <= options.min_stepsize*1.001_num) { // cannot achieve better accuracy
            step_success= false;
            break;
        }
        step_size = std::max(dr/4._num, options.min_stepsize); // otherwise try again with smaller stepsize
    }
    current_step.first = r;
    current_step.second = y;

    return step_success;
}


/* This function integrates an ODE system according to the RKF45 algorithm
 *
 * The system is described by dy_dr, the integration starts at r0 and tries to reach r_end
 * The initial values are given by y0. The params pointer can be arbitrary and is passed on to dy_dr
 * The initial and last step are saved in results (unless options.save_intermediate == true)
 * and any number of events can be tracked throughout the evolution with the events vector */
int integrator::RKF45(ODE_system dy_dr, const NUMERIC r0, const vector y0, const NUMERIC r_end, const void* params,
                            std::vector<step>& results, std::vector<Event>& events, const IntegrationOptions& options)
{
    NUMERIC step_size = options.max_stepsize;        // initial stepsize
    vector dy;
    integrator::step current_step = std::make_pair(r0, y0);

    // clear passed arguments
    results.clear();
    if(options.clean_events) {
        for(auto it = events.begin(); it != events.end(); ++it)
            it->reset();
    }

    // reserve space if necessary
    if(options.save_intermediate) {
        results.reserve(options.max_step);
        results.push_back(current_step);
    }

    if(options.verbose > 0)
        std::cout << "starting integration at r=" << current_step.first << std::endl;

    // begin integration
    int i = 0, stop = 0;
    while(true) {

        //bool step_success = RKF45_step(dy_dr, current_step.first, step_size, current_step.second, params, options);
        bool step_success = RKF45_step_event_tester(dy_dr, current_step, step_size, params, events,  options);
        dy = dy_dr(current_step.first, current_step.second, params);
        if(options.verbose > 2)
            std::cout  << "rkf45 step: r=" <<  current_step.first << ", dr= " << step_size << ", y=" << current_step.second << ", dy=" << dy <<  std::endl;

        if(options.save_intermediate)
            results.push_back(current_step);

        // iterate through all defined events and check them (happens every iteration step)
        for(auto it = events.begin(); it != events.end(); ++it) {
            if(it->condition(current_step.first, step_size, current_step.second, dy, params)) {
                if(!it->active) {   // these events only trigger when the condition wasn't previously active
                    it->active = true;
                    it->steps.push_back(current_step);  // add the current values of the ODE vars to the event object
                    if(it->stopping_condition) {     // if the event is defined as a stopping condition, we stop the iteration (see below)
                        stop = event_stopping_condition;
                        if(options.verbose > 0)
                            std::cout << "stopping condition: " << it->name << " triggered - ";
                    }
                }
            } else
                it->active = false;
        }

        // check additional stopping conditions
        if(current_step.first > r_end)     // max r reached
            stop = endpoint_reached;
        if(i > options.max_step)            // max number of steps reached
            stop = iteration_number_exceeded;
        if(!step_success)                   // problem in the last step (e.g. stepsize too small and/or truncation error too large)
            stop = stepsize_underflow;
        if(stop) {
            if(!options.save_intermediate)  // save the last step if not already done
                results.push_back(current_step);
            if(options.verbose > 0)
                std::cout << "stopped integration with code " << stop << " at r=" << current_step.first << std::endl;
            return stop;    // exit the integrator if a stopping condition was reached
        }
        i++;
    }
}

/* cumulative trapezoid integration for the pair x,y
 * output is saved in res */
void integrator::cumtrapz(const std::vector<NUMERIC>& x, const std::vector<NUMERIC>& y, std::vector<NUMERIC>& res) {
    if( x.size() != y.size() || x.size() == 0)
        return;
    res = std::vector<NUMERIC>(x.size(), 0._num);

    for(unsigned int i = 1; i < x.size(); i++) {
        res[i] = (x[i]-x[i-1]) * (y[i] + y[i-1])/2._num;
        res[i] += res[i-1];
    }
}

std::ostream& integrator::operator <<(std::ostream& o, const integrator::step& s) {
    return o << "r = " << s.first << ", y = " << s.second;
}
