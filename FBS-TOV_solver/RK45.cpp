#include "RK45.hpp"

// Runge-Kutta-Fehlberg stepper for one integration step:
// for the scheme in detail, see https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method -> COEFFICIENTS FOR RK4(5), FORMULA 2 Table III in Fehlberg
void RK45::RK45_Fehlberg(ODE_system dy_dt, double &r, double &dr , vec5d& V, const void *params)
{
    // variables for error estimation:
    const double target_error = 1e-9;   // what truncation error should we allow/aim for?
    double max_stepsize = 1e-3;   // maximum stepsize in units of dr
    const double min_stepsize = 1e-5;   // minimum stepsize in untis of dr (to prevent deadlock due to decreasing stepsize)

    // intermediate steps:
    vec5d k1, k2, k3, k4, k5, k6;

    // solve until we find an acceptable stepsize/error:
    while (true) {

        // runge kutta Fehlberg steps for all ODEs:
        k1 = dr * dy_dt(r, V, params);
        k2 = dr * dy_dt(r + 1.0 / 4.0 * dr, V + 1.0 / 4.0 * k1, params);
        k3 = dr * dy_dt(r + 3.0 / 8.0 * dr, V + 3.0 / 32.0 * k1 + 9.0 / 32.0 * k2, params);
        k4 = dr * dy_dt(r + 12.0 / 13.0 * dr, V + 1932.0 / 2197.0 * k1 - 7200.0 / 2197.0 * k2 + 7296.0 / 2197.0 * k3, params);
        k5 = dr * dy_dt(r + dr, V + 439.0 / 216.0 * k1 - 8.0 * k2 + 3680.0 / 513.0 * k3 - 845.0 / 4104.0 * k4, params);
        k6 = dr * dy_dt(r + 1.0 / 2.0 * dr, V - 8.0 / 27.0 * k1 + 2.0 * k2 - 3544.0 / 2565.0 * k3 + 1859.0 / 4104.0 * k4 - 11.0 / 40.0 * k5, params);

        // 4th and 5th order accurate steps:
		vec5d dV_O4 = 25.0 / 216.0 * k1 + 1408.0 / 2565.0 * k3 + 2197.0 / 4104.0 * k4 - 1.0 / 5.0 * k5; // + O(x^5)
        vec5d dV_O5 = 16.0 / 135.0 * k1 + 6656.0 / 12825.0 * k3 + 28561.0 / 56430.0 * k4 - 9.0 / 50.0 * k5 + 2.0 / 55.0 * k6; // + O(x^6)

		// approximating the truncation error:
        double truncation_error = vec5d::max(dV_O5 - dV_O4) * dr; // inf-Norm
		//double truncation_error = vec5d::dot(dV_O5 - dV_O4, dV_O5 - dV_O4) * dr; // 2-Norm

		// Updating the stepsize:

        if( std::isnan(dV_O5[0]) || std::isnan(dV_O5[1]) || std::isnan(dV_O5[2]) || std::isnan(dV_O5[3]) || std::isnan(dV_O5[4])) {
        std::cout << "Nan found" << std::endl;
                exit(1);}

        //if (V[4] < 1e-15) max_stepsize = 0.1;   // allow a larger stepsize outside of the star (defined by P=0) (mainly we want an accurate result for the star radius)

		if (truncation_error > target_error) {

            if (dr < min_stepsize) {
                // error is not acceptable but the stepsize cannot get any smaller:
                // write new values for r and the evolved quantities
                r += dr;
                V += dV_O5;

                dr = min_stepsize; // ensure that the stepsize never gets too small
                std::cout << "Error in RK45_Fehlberg(): Minimal stepsize underflow at r = " << r << ", dr = "<< dr << std::endl;  // error message for debugging
                std::cout << "Truncation error = "<< truncation_error << std::endl;
                std::cout << V[0] << " " << V[1] << " " << V[2] << " " << V[3] << " " << V[4] << " " << V[5] << std::endl;
                return;
            }

			dr *= 0.5;     // truncation error is too large. We repeat the iteration with smaller stepsize
			continue;
		}
		else {
			// error is small enough and therefore we can use this step. To save computation time, we then increase the stepsize
            r += dr;
            V += dV_O5;

            dr *= 2.0;
			if (dr > max_stepsize) dr = max_stepsize;   // enforce maximal stepsize

			return;
		}

	    return;  // stop this iteration step if the required accuracy and step size was acheived
	}
}
