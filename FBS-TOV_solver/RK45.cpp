#include "RK45.hpp"

namespace ublas = boost::numeric::ublas;


bool integrator::RKF45_step(ODE_system dy_dt, double &r, double &dr, vector& y, const void* params, const IntegrationOptions& options)
{
    // intermediate steps:
    int n = y.size();
    vector k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), dV_04(n), dV_05(n);

    // solve until we find an acceptable stepsize/error:
    while (true) {

        // runge kutta Fehlberg steps for all ODEs:
        k1 = dr * dy_dt(r, y, params);
        k2 = dr * dy_dt(r + 1.0 / 4.0 * dr, y + 1.0 / 4.0 * k1, params);
        k3 = dr * dy_dt(r + 3.0 / 8.0 * dr, y + 3.0 / 32.0 * k1 + 9.0 / 32.0 * k2, params);
        k4 = dr * dy_dt(r + 12.0 / 13.0 * dr, y + 1932.0 / 2197.0 * k1 - 7200.0 / 2197.0 * k2 + 7296.0 / 2197.0 * k3, params);
        k5 = dr * dy_dt(r + dr, y + 439.0 / 216.0 * k1 - 8.0 * k2 + 3680.0 / 513.0 * k3 - 845.0 / 4104.0 * k4, params);
        k6 = dr * dy_dt(r + 1.0 / 2.0 * dr, y - 8.0 / 27.0 * k1 + 2.0 * k2 - 3544.0 / 2565.0 * k3 + 1859.0 / 4104.0 * k4 - 11.0 / 40.0 * k5, params);

        // 4th and 5th order accurate steps:
		dV_04 = 25.0 / 216.0 * k1 + 1408.0 / 2565.0 * k3 + 2197.0 / 4104.0 * k4 - 1.0 / 5.0 * k5; // + O(x^5)
        dV_05 = 16.0 / 135.0 * k1 + 6656.0 / 12825.0 * k3 + 28561.0 / 56430.0 * k4 - 9.0 / 50.0 * k5 + 2.0 / 55.0 * k6; // + O(x^6)

        if( vector::is_nan(dV_05) || vector::is_nan(dV_04)) {
        std::cout << "Nan found" << dV_05 << dV_04 << k1 << k2 << k3 << k4 << k5 << k6 << "\n  for r=" << r << ", dr=" << dr <<  std::endl;
                exit(1);}

        // approximating the truncation error:
        double truncation_error = ublas::norm_inf(dV_05 - dV_04) * dr; // inf-Norm
		//double truncation_error = ublas::norm_2(dV_05 - dV_04) * dr; // 2-Norm


        //if (y[4] < 1e-15) max_stepsize = 0.1;   // allow a larger stepsize outside of the star (defined by P=0) (mainly we want an accurate result for the star radius)

		if (truncation_error > options.target_error) {
			dr *= 0.5;     // truncation error is too large. We repeat the iteration with smaller stepsize

            if (dr < options.min_stepsize) {
                // error is not acceptable but the stepsize cannot get any smaller:
                // write new values for r and the evolved quantities
                r += dr;
                y += dV_05;

                dr = options.min_stepsize; // ensure that the stepsize never gets too small
                std::cout << "Error in RKF45_step(): Minimal stepsize underflow at r = " << r << ", dr = "<< dr << std::endl;  // error message for debugging
                std::cout << "Truncation error = "<< truncation_error << std::endl;
                std::cout << y << std::endl;
                return false;
            }

			continue;
		}
		else {
			// error is small enough and therefore we can use this step. To save computation time, we then increase the stepsize
            r += dr;
            y += dV_05;

            dr *= 2.0;
			if (dr > options.max_stepsize) dr = options.max_stepsize;   // enforce maximal stepsize

			return true;
		}

	    return true;  // stop this iteration step if the required accuracy and step size was acheived
	}
}

int integrator::RKF45(ODE_system dy_dt, const double r0, const vector y0, const double r_end, const void* params,
                            std::vector<step>& results, std::vector<Event>& events, const IntegrationOptions& options)
{
    double step_size = options.max_stepsize;        // initial stepsize
    vector dy;
    integrator::step current_step = std::make_pair(r0, y0);

    // clear passed arguments
    results.clear();
    for(auto it = events.begin(); it != events.end(); ++it)
        it->reset();

    // reserve space if necessary
    if(options.save_intermediate) {
        results.reserve(options.max_step);
        results.push_back(current_step);
    }

    // begin integration
    int i = 0, stop = 0;
    while(true) {

        bool step_success = RKF45_step(dy_dt, current_step.first, step_size, current_step.second, params, options);
        //std::cout << "rkf45 step: r=" <<  current_step.first << ", dr= " << step_size << ", y=" << current_step.second << std::endl;
        dy = dy_dt(current_step.first, current_step.second, params);

        if(options.save_intermediate)
            results.push_back(current_step);

        for(auto it = events.begin(); it != events.end(); ++it) {
            if(it->condition(current_step.first, step_size, current_step.second, dy, params)) {
                if(!it->active) {
                    it->active = true;
                    it->steps.push_back(current_step);
                    if(it->stopping_condition)
                        stop = event_stopping_condition;
                }
            } else
                it->active = false;
        }

        if (current_step.first > r_end)
            stop = endpoint_reached;
        if(i > options.max_step)
            stop = iteration_number_exceeded;
        if(!step_success)
            stop = stepsize_underflow;
        if(stop) {
            if(!options.save_intermediate)
                results.push_back(current_step);
            std::cout << "stopped integration with code " << stop << " at r=" << current_step.first << std::endl;
            return stop;
        }
        i++;
    }
}
