#include "RK45.hpp"

namespace ublas = boost::numeric::ublas;


void integrator::RKF45_step(ODE_system dy_dt, double &r, double &dr, vector& y, const void* params, const double target_error, const double max_stepsize, const double min_stepsize)
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

		// approximating the truncation error:
        double truncation_error = ublas::norm_inf(dV_05 - dV_04) * dr; // inf-Norm
		//double truncation_error = vector::dot(dV_05 - dV_04, dV_05 - dV_04) * dr; // 2-Norm

		// Updating the stepsize:

        if( std::isnan(dV_05[0]) || std::isnan(dV_05[1]) || std::isnan(dV_05[2]) || std::isnan(dV_05[3]) || std::isnan(dV_05[4])) {
        std::cout << "Nan found" << std::endl;
                exit(1);}

        //if (y[4] < 1e-15) max_stepsize = 0.1;   // allow a larger stepsize outside of the star (defined by P=0) (mainly we want an accurate result for the star radius)

		if (truncation_error > target_error) {

            if (dr < min_stepsize) {
                // error is not acceptable but the stepsize cannot get any smaller:
                // write new values for r and the evolved quantities
                r += dr;
                y += dV_05;

                dr = min_stepsize; // ensure that the stepsize never gets too small
                std::cout << "Error in RKF45_step(): Minimal stepsize underflow at r = " << r << ", dr = "<< dr << std::endl;  // error message for debugging
                std::cout << "Truncation error = "<< truncation_error << std::endl;
                std::cout << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << std::endl;
                return;
            }

			dr *= 0.5;     // truncation error is too large. We repeat the iteration with smaller stepsize
			continue;
		}
		else {
			// error is small enough and therefore we can use this step. To save computation time, we then increase the stepsize
            r += dr;
            y += dV_05;

            dr *= 2.0;
			if (dr > max_stepsize) dr = max_stepsize;   // enforce maximal stepsize

			return;
		}

	    return;  // stop this iteration step if the required accuracy and step size was acheived
	}
}

int integrator::RKF45(ODE_system dy_dt, const double r0, const vector y0, const double r_end, const void* params, const int max_step,
                            std::vector<step>& results, std::vector<Event>& events, const bool save_intermediate)
{
    double step_size = 1e-5;        // initial stepsize

    integrator::step current_step = std::make_pair(r0, y0);
    results.clear();
    if(save_intermediate) {
        results.reserve(max_step);
        results.push_back(current_step);
    }

    int i = 0, stop = 0;
    while(true) {

        RKF45_step(dy_dt, current_step.first, step_size, current_step.second, params);

        if(save_intermediate)
            results.push_back(current_step);

        for(auto it = events.begin(); it != events.end(); ++it) {
            if(it->condition(current_step.first, current_step.second, params)) {
                if(!it->active) {
                    it->active = true;
                    it->steps.push_back(current_step);
                    if(it->stopping_condition)
                        stop = 3;
                }
            } else
                it->active = false;
        }

        if (current_step.first > r_end)
            stop = 2;
        if(i > max_step)
            stop = 1;
        if(stop) {
            if(!save_intermediate)
                results.push_back(current_step);
            std::cout << "stopped integration with code " << stop << " at r=" << current_step.first << std::endl;
            return stop;
        }
        i++;
    }
}
