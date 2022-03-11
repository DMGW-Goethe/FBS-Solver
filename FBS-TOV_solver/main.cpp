#include <iostream> // for i/o e.g. std::cout etc.
#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <fstream>	// file streams
// #include <mpi.h>  // to use MPI (for parallelization) -> needs special compile rules!

#include "vec5d.hpp"    // include custom 5-vector class
#include "eostable.hpp" // include eos container class

// --------------------------------------------------------------------

// polytropic EOS: p=kappa*rho**Gamma
void EOS_polytrope(double& myrho, double& myepsilon, const double P) {

    const double kappa = 100.; // polytropic constant
    const double Gamma = 2.;  // adiabatic index

    // for these hydrodynamic relations see e.g. "Relativistic hydrodynamics, Rezzolla and Zanotti" chapter 2.4.7 pages 118+119
    myrho = std::pow(P / kappa, 1./Gamma); // update restmass density according to polytropic EOS
    myepsilon = kappa*std::pow(myrho, Gamma - 1.) / (Gamma - 1.);
}


// define the system of equations:
// implemented as in https://arxiv.org/abs/2110.11997 eq. 7-11
// parameters:  radius r (evolution parameter), vec5d containing all evolved vars, mass mu, self-interaction parameter lambda, frequency omega:  
vec5d ODE_system(const double r, const vec5d& vars, eostable& myEOS, const double mu, const double lambda, const double omega) {

	// rename input variables for simpler use:
	const double a = vars[0]; const double alpha = vars[1]; const double Phi = vars[2]; const double Psi = vars[3]; double P = vars[4];

    // define hydrodynamic quantities:
    double rho = 1.;      // restmass density, must be set using EOS
    double epsilon = 1.;  // specific energy denstiy, must be set either through EOS or hydrodynamic relations
                          // epsilon is related to the total energy density "e" by: e = rho*(1+epsilon)
    
    if(P < 0.) P = 0.;  // need this to prevent NaN errors...

    // apply the EOS:
    myEOS.callEOS(rho, epsilon, P); // change rho and epsilon by pointer using EOS member function

    if( std::isnan(vars[0]) || std::isnan(vars[1]) || std::isnan(vars[2]) || std::isnan(vars[3]) || std::isnan(vars[4])) {
        std::cout << "Nan found" << std::endl;
                exit(1);} 

	// compute the ODEs:
	double da_dr = 0.5* a * ( (1.-a*a) / r + 4.*M_PI*r*( (omega*omega/ alpha*alpha + mu*mu + 0.5*lambda*Phi*Phi )*a*a*Phi*Phi + Psi*Psi + 2.*a*a*rho*(1.+epsilon) ) ); // a (g_rr component) ODE
	double dalpha_dr = 0.5* alpha * ( (a*a-1.) / r + 4.*M_PI*r*( (omega*omega/ alpha*alpha - mu*mu - 0.5*lambda*Phi*Phi )*a*a*Phi*Phi + Psi*Psi + 2.*a*a*P ) ); // alpha (g_tt component) ODE
	double dPhi_dr = Psi; // Phi ODE
	double dPsi_dr = -( 1. + a*a - 4.*M_PI*r*r*a*a*( mu*mu*Phi*Phi + 0.5*lambda*Phi*Phi*Phi*Phi + rho*(1.+epsilon) - P ))*Psi/r - (omega*omega/ alpha*alpha - mu*mu - lambda*Phi*Phi )*a*a*Phi*Phi; // Psi ODE
	double dPres_dr = -(rho*(1.+epsilon) + P)*dalpha_dr/alpha; // Pressure ODE

	// write the ODE values into output vector:
	return vec5d(da_dr, dalpha_dr, dPhi_dr, dPsi_dr, dPres_dr);
}


// Runge-Kutta-Fehlberg stepper for one integration step:
// for the scheme in detail, see https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method -> COEFFICIENTS FOR RK4(5), FORMULA 2 Table III in Fehlberg
void RK45_Fehlberg(double &r, double &dr , vec5d& V, eostable& myEOS, const double mu, const double lambda, const double omega)
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
        k1 = dr * ODE_system(r, V, myEOS, mu, lambda, omega);
        k2 = dr * ODE_system(r + 1.0 / 4.0 * dr, V + 1.0 / 4.0 * k1, myEOS, mu, lambda, omega);
        k3 = dr * ODE_system(r + 3.0 / 8.0 * dr, V + 3.0 / 32.0 * k1 + 9.0 / 32.0 * k2, myEOS, mu, lambda, omega);
        k4 = dr * ODE_system(r + 12.0 / 13.0 * dr, V + 1932.0 / 2197.0 * k1 - 7200.0 / 2197.0 * k2 + 7296.0 / 2197.0 * k3, myEOS, mu, lambda, omega);
        k5 = dr * ODE_system(r + dr, V + 439.0 / 216.0 * k1 - 8.0 * k2 + 3680.0 / 513.0 * k3 - 845.0 / 4104.0 * k4, myEOS, mu, lambda, omega);
        k6 = dr * ODE_system(r + 1.0 / 2.0 * dr, V - 8.0 / 27.0 * k1 + 2.0 * k2 - 3544.0 / 2565.0 * k3 + 1859.0 / 4104.0 * k4 - 11.0 / 40.0 * k5, myEOS, mu, lambda, omega);
		
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


// integrating the system of equations using the shooting method:
// for shoooting method see e.g. https://sdsawtelle.github.io/blog/output/shooting-method-ode-numerical-solution.html
// 2D array to store all physical values, initial values, mu, lambda
// returns all physical values in a 2D array, also R and M of the star
template <typename Matrix2D>
void Shooting_integrator_save_intermediate(Matrix2D& myarray, const unsigned Nmax_iterations, double& Rfermionic, double& Mtot, vec5d& init_vars, eostable& myEOS, const double mu, const double lambda) {

    // define the evolution variables
    double omega_bar = 1.0;             // initial guess for omega and running variable
    double omega_0bestguess = omega_bar; // this will be the best guess for omega
    double delta_omega = 1e-5;           // to track the quality of the derivative for finite-difference
    const double root_desired_error = 1e-4;    // desired accuracy for the root F(omega)=0

    const double init_r = 1e-10;        // initial radius (cannot be zero due to coordinate singularity)
    const double init_dr = 1e-5;        // initial stepsize

    // declare the variables to be able to use them in the whole function scope
    double my_r = init_r;
    double my_dr = init_dr;
    vec5d ODE_vars = init_vars;
    vec5d ODE_vars2 = init_vars;

    // fill the array with the init values:
    myarray[0][0] = my_r;           // radius
	myarray[0][1] = ODE_vars[0];    // a (g_rr component)
    myarray[0][2] = ODE_vars[1];    // alpha (g_tt component)
    myarray[0][3] = ODE_vars[2];    // Phi (scalar field)
	myarray[0][4] = ODE_vars[3];    // Psi (derivative of Phi)
	myarray[0][5] = ODE_vars[4];    // P (pressure)

    // integrate one TOV star until a stable configuration has been reached (including shooting method):
    // outer shooting method loop:
    while(true) {

        // reset the initial values before each ODE itegration (only omega should change each iteration):
        ODE_vars = init_vars;
        ODE_vars2 = init_vars; 
        my_r = init_r;
        my_dr = init_dr;
        Rfermionic = 0.0;
        Mtot = 0.0;

        // inner ODE integration loop:
        for (unsigned i = 1; i < Nmax_iterations+1; ++i) {
            
            // step one ODE iteration:
            RK45_Fehlberg(my_r, my_dr, ODE_vars, myEOS, mu, lambda, omega_bar);
            
            // save the results of this step into a suitable 2D array/container:
            myarray[i][0] = my_r;           // radius
			myarray[i][1] = ODE_vars[0];    // a (g_rr component)
			myarray[i][2] = ODE_vars[1];    // alpha (g_tt component)
            myarray[i][3] = ODE_vars[2];    // Phi (scalar field)
			myarray[i][4] = ODE_vars[3];    // Psi (derivative of Phi)
			myarray[i][5] = ODE_vars[4];    // P (pressure)

            // check for stopping condition (Abbruchbedingungen):
            if ( (ODE_vars[4] < 1e-15) && (Rfermionic < 0.01) ) { // pressure is zero -> fermionic radius was found
                Rfermionic = my_r;
            }
            // radius is large (so that we include the whole bosonic scalar field as well):
            if (my_r > 100.) {
                
                // fill all remaining array elements with the last integrated value:
                for (unsigned j = i+1; j<Nmax_iterations+1; ++j) {
				    myarray[j][0] = my_r;           // radius
				    myarray[j][1] = ODE_vars[0];    // a (g_rr component)
				    myarray[j][2] = ODE_vars[1];    // alpha (g_tt component)
                    myarray[j][3] = ODE_vars[2];    // Psi (derivative of Phi)
				    myarray[j][4] = ODE_vars[3];    // Phi (scalar field)
				    myarray[j][5] = ODE_vars[4];    // P (pressure)
			    }
                break;} // radius is 100M
        }
        // reset all init values before integrating the second time (to perform the finite difference step for omega):
        my_r = init_r;
        my_dr = init_dr;
        // inner ODE integration loop second solution (finite difference in omega was added):
        for (unsigned i = 1; i < Nmax_iterations+1; ++i) {
            // step two ODE iteration:
            RK45_Fehlberg(my_r, my_dr, ODE_vars2, myEOS, mu, lambda, omega_bar+delta_omega); // for the NR method relative to omega:
            
            // check for stopping condition (Abbruchbedingungen):
            // integrate until radius is large (so that we include the whole bosonic scalar field as well):
            if (my_r > 100.) {break;} // radius is 100M
        }
        
        // take result from the integrator and try to obtain a new value for omega:
        // Newton-Raphson method to step in the right direction for omega
        omega_0bestguess =  omega_bar - ODE_vars[2] * delta_omega / ( ODE_vars2[2] - ODE_vars[2] );

        std::cout << omega_bar << std::endl;

		if ( std::abs( ODE_vars2[2] ) < root_desired_error) { // stop the loop once a satisfactory level of accurace has been found
			
            // exit the shooting integrator as the wanted accuracy has been found
            break;
		}
		else {
			omega_bar = omega_0bestguess; // let x_bar become the starting point for the next iteration step
		}

    }

    // optional: calculate the last values here...
    // e.g. calculate star mass and fermionic radius
    // Rfermionic was already computed (see above)
    // compute total mass (bosonic+fermionic) as " Mtot = lim_{r->inf} r/2 * (1 - 1/g_rr) ", see https://arxiv.org/abs/2110.11997 eq. 13
    Mtot = 0.5 * my_r * ( 1.0 - 1.0/(ODE_vars[0]*ODE_vars[0]) );
}


// integrating the system of equations using the shooting method:
// for shoooting method see e.g. https://sdsawtelle.github.io/blog/output/shooting-method-ode-numerical-solution.html
// parameters: star Radius, Mass, initial values, mu, lambda
// returns only the Mass and fermionic radius of ONE star
void Shooting_integrator_nosave_intermediate(double& Rfermionic, double& Mtot, vec5d& init_vars, eostable& myEOS, const double mu, const double lambda) {

    // define the evolution variables
    double omega_bar = 1.0;             // initial guess for omega and running variable
    double omega_0bestguess = omega_bar; // this will be the best guess for omega
    double delta_omega = 1e-5;           // to track the quality of the derivative for finite-difference
    const double root_desired_error = 1e-4;    // desired accuracy for the root F(omega)=0

    const double init_r = 1e-10;        // initial radius (cannot be zero due to coordinate singularity)
    const double init_dr = 1e-5;        // initial stepsize

    // declare the variables to be able to use them in the whole function scope
    double my_r = init_r;
    double my_dr = init_dr;
    vec5d ODE_vars = init_vars;
    vec5d ODE_vars2 = init_vars;
    
    // integrate one TOV star until a stable configuration has been reached (including shooting method):
    // outer shooting method loop:
    while(true) {

        // reset the initial values before each ODE itegration (only omega should change each iteration):
        ODE_vars = init_vars;
        ODE_vars2 = init_vars; 
        my_r = init_r;
        my_dr = init_dr;
        Rfermionic = 0.0;
        Mtot = 0.0;
        // inner ODE integration loop first solution:
        while (true) {
            // step one ODE iteration:
            RK45_Fehlberg(my_r, my_dr, ODE_vars, myEOS, mu, lambda, omega_bar);
            
            // detect the fermionic radius defined by the point where P=0.
            if ( (ODE_vars[4] < 1e-14) && (Rfermionic < 0.01) ) { // pressure is zero -> fermionic radius was found
                Rfermionic = my_r;
            }
            // check for stopping condition (Abbruchbedingungen):
            // integrate until radius is large (so that we include the whole bosonic scalar field as well):
            if (my_r > 100.) {break;} // radius is 100M
            //std::cout << "one_Fehlberg_iteration done loop1! at r=" << my_r << " fermi rad:" << Rfermionic <<std::endl;
            //std::cout << "(a,alpha,Phi,Psi,P)=" << ODE_vars[0] <<" "<< ODE_vars[1] <<" "<< ODE_vars[2] <<" "<< ODE_vars[3] <<" "<< ODE_vars[4] << std::endl;
           // if( std::isnan(ODE_vars[0])) {exit(1);}    // nan error
        }
        // reset all init values before integrating the second time (to perform the finite difference step for omega):
        my_r = init_r;
        my_dr = init_dr;
        // inner ODE integration loop second solution (finite difference in omega was added):
        while (true) {
            // step two ODE iteration:
            RK45_Fehlberg(my_r, my_dr, ODE_vars2, myEOS, mu, lambda, omega_bar+delta_omega); // for the NR method relative to omega
            
            // check for stopping condition (Abbruchbedingungen):
            // integrate until radius is large (so that we include the whole bosonic scalar field as well):
            if (my_r > 100.) {break;} // radius is 100M
        }

    //if (std::abs( ODE_vars2[2]) < root_desired_error) { break;}

        // take result from the integrator and try to obtain a new value for omega:
        // Newton-Raphson method to step in the right direction for omega
        omega_0bestguess =  omega_bar - ODE_vars[2] * delta_omega / ( ODE_vars2[2] - ODE_vars[2] );

		if ( std::abs( ODE_vars2[2] ) < root_desired_error) { // stop the loop once a satisfactory level of accuracy has been found
			
            // exit the shooting integrator as the wanted accuracy has been found
            break;
		}
		else {
			omega_bar = omega_0bestguess; // let the best guess become the starting point for the next iteration step
		}

        std::cout << "one_omega_iteration done! omega = " << omega_0bestguess << std::endl;
    }

    // optional: calculate the last values here...
    // e.g. calculate star mass and fermionic radius
    // Rfermionic was already computed (see above)
    // compute total mass (bosonic+fermionic) as " Mtot = lim_{r->inf} r/2 * (1 - 1/g_rr) ", see https://arxiv.org/abs/2110.11997 eq. 13
    //std::cout << "compute Mtot " << ODE_vars[0] <<std::endl;
    Mtot = 0.5 * my_r * ( 1.0 - 1.0/(ODE_vars[0]*ODE_vars[0]) );
}


// saves the calculated ODEs paramaters in a text file
template <typename Mat2D>
void save_data_ode(Mat2D& array, const unsigned length, std::string filename) {
	
	std::ofstream img;
	img.open(filename);

	if(img.is_open()) {

		img << "# r\t     a\t    alpha\t    Phi\t    Psi\t    P" << std::endl;

		for(unsigned i = 0; i < length; ++i) {

			img <<std::fixed<<std::setprecision(10)<< array[i][0] << " " << array[i][1] << " " << array[i][2] << " " << array[i][3] << " " << array[i][4] << " " << array[i][5] << std::endl;
		}
	}
	img.close();
}


// saves the MR-values in a text file
template <typename Mat2D>
void save_data_MR(Mat2D& array, const unsigned length, std::string filename) {
	
	std::ofstream img;
	img.open(filename);

	if(img.is_open()) {

		img << "# R [km]\t M [M_sun]\t rho_c [code units]" << std::endl;

		for(unsigned i = 0; i < length; ++i) {

			img <<std::fixed<<std::setprecision(10)<< array[i][0]*1.476625061 << " " << array[i][1] << " " << array[i][2] << std::endl;
		}
	}
	img.close();
}

int main() {

    // for ONE star:
    const unsigned N_steps = 100000;   // number of iteration steps in the ODE solver (needs to be fie-tuned later). Only use if we want to save all intermediate values!
	double ODEdata[N_steps+1][6] = {0.};  // container for the data in the ODE integrator. 

    // for the MR diagram (multiple stars):
    const unsigned Nstars = 10;
    double MRarray[Nstars][3] = {0.};

    // define some global values:
    double mu = 0.0;        // DM mass
    double lambda = 0.0;    //self-interaction parameter
    double R_fermi = 0.0; double M_total = 0.0; // fermionic radius of star and total mass of star

    double rho_c = 2e-4; // central density

    // integrate ONE star and save all intermediate values into a txt file:
    // a, alpha, Phi, Psi, P(rho)
    //vec5d inits(1.0, 1.0, 1e-20, 1e-20, 100*std::pow(10*rho_c, 2.) );
    // solve the ODEs WITH saing all intermediate steps.
    //Shooting_integrator_save_intermediate(ODEdata, N_steps, R_fermi, M_total, inits, mu, lambda);
    //save_data_ode(ODEdata, N_steps, "ODE_data_full.txt");
    // print the results for now (later save them into a txt file).
    //std::cout << "Star with rho_c = " << rho_c << ": radius = " << R_fermi << " [M], mass = " << M_total << " [M_sun]" << std::endl; 

    // try the new tabulated EOS system:
    eostable myDD2("DD2_eos.table");
    eostable myPolytrope; // default constructor => polytrope

    //vec5d inits(1.0, 1.0, 1e-20, 1e-20, 100*std::pow(10*rho_c, 2.) );
    //Shooting_integrator_nosave_intermediate(R_fermi, M_total, inits, myDD2, mu, lambda);
    

    // uncomment this if you want to compute a full MR-diagram

    std::cout << "start loop" << std::endl;

    // create a simple MR-diagram:
    for (unsigned i = 0; i < Nstars; ++i) {

        double rho_start = (i+1)*rho_c;

        // a, alpha, Phi, Psi, P(rho)
        vec5d inits(1.0, 1.0, 1e-20, 1e-20, myPolytrope.get_P_from_rho(rho_start) );
        Shooting_integrator_nosave_intermediate(R_fermi, M_total, inits, myPolytrope, mu, lambda);
        //Shooting_integrator_nosave_intermediate(R_fermi, M_total, inits, myDD2, mu, lambda);

        MRarray[i][0] = R_fermi;
        MRarray[i][1] = M_total;
        MRarray[i][2] = rho_start;

       std::cout << i << std::endl;
    }

    std::cout << "end loop" << std::endl;

    // print the results:
    for (unsigned j = 0; j < Nstars; ++j) {
        //std::cout << j << std::endl;
        std::cout << "R = " << MRarray[j][0]*1.476625061 << " [km], M = " << MRarray[j][1] << " [M_sun]" << std::endl;
    }

    // save MR data in text file:
    save_data_MR(MRarray, Nstars, "MRtest_DD2eos.txt");


    // print the results for now (later save them into a txt file).
    //std::cout << "Star with rho_c = " << rho_c << ": radius = " << R_fermi << " [M], mass = " << M_total << " [M_sun?]" << std::endl; 


    return 0;
}
