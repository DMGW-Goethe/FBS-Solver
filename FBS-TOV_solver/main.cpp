#include <iostream> // for i/o e.g. std::cout etc.
#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <fstream>	// file streams
#include <memory>
// #include <mpi.h>  // to use MPI (for parallelization) -> needs special compile rules!

#include "vector.hpp"    // include custom 5-vector class
#include "RK45.hpp"
#include "eos.hpp" // include eos container class

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
// parameters:  radius r (evolution parameter), vector containing all evolved vars, mass mu, self-interaction parameter lambda, frequency omega:
vector dy_dt(const double r, const vector& vars, EquationOfState& myEOS, const double mu, const double lambda, const double omega) {

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
	return vector({da_dr, dalpha_dr, dPhi_dr, dPsi_dr, dPres_dr});
}

struct SystemParameters {
    std::shared_ptr<EquationOfState> EOS;
    double mu;
    double lambda;
    double omega;
};

vector ODE_system(const double r, const vector& y, const void* params) {

    SystemParameters* sp = (SystemParameters*)params;
    return dy_dt(r, y, *(sp->EOS), sp->mu, sp->lambda, sp->omega);
}


// integrating the system of equations using the shooting method:
// for shoooting method see e.g. https://sdsawtelle.github.io/blog/output/shooting-method-ode-numerical-solution.html
// parameters: star Radius, Mass, initial values, mu, lambda
// returns only the Mass and fermionic radius of ONE star
void Shooting_integrator_nosave_intermediate(double& Rfermionic, double& Mtot, vector& init_vars, SystemParameters* sp) {

    // define the evolution variables
    double omega_bar = 1.0;             // initial guess for omega and running variable
    double omega_0bestguess = omega_bar; // this will be the best guess for omega
    double delta_omega = 1e-5;           // to track the quality of the derivative for finite-difference
    const double root_desired_error = 1e-4;    // desired accuracy for the root F(omega)=0

    const double init_r = 1e-10;        // initial radius (cannot be zero due to coordinate singularity)
    const double init_dr = 1e-5;        // initial stepsize
    const double r_end = 100.;
    int max_step = 1000000;

    // declare the variables to be able to use them in the whole function scope
    double my_r = init_r;
    double my_dr = init_dr;
    vector ODE_vars = init_vars;
    vector ODE_vars2 = init_vars;

    auto Rf_condition = [](const double r, const vector y, const void*params) { return y[4] < 1e-14; };
    std::vector<RK45::step> results, events;

    // integrate one TOV star until a stable configuration has been reached (including shooting method):
    // outer shooting method loop:
    while(true) {

        // reset the initial values before each ODE itegration (only omega should change each iteration):
        Rfermionic = 0.0;
        Mtot = 0.0;
        // inner ODE integration loop first solution:
        sp->omega = omega_bar;
        int res =  RK45::RKF45(&ODE_system, init_r, init_vars, r_end, (void*) sp, max_step,  results,  events, false, Rf_condition);
        ODE_vars = results[0].second; my_r = results[0].first;

        Rfermionic = events.at(0).first;
        // second solution
        sp->omega = omega_bar + delta_omega;
        res =  RK45::RKF45(&ODE_system, init_r, init_vars, r_end, (void*) sp, max_step,  results,  events);
        ODE_vars2 = results[0].second;

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
    SystemParameters sp;
    sp.mu = 0.0;        // DM mass
    sp.lambda = 0.0;    //self-interaction parameter
    double R_fermi = 0.0; double M_total = 0.0; // fermionic radius of star and total mass of star

    double rho_c = 2e-4; // central density

    // integrate ONE star and save all intermediate values into a txt file:
    // a, alpha, Phi, Psi, P(rho)
    //vector inits(1.0, 1.0, 1e-20, 1e-20, 100*std::pow(10*rho_c, 2.) );
    // solve the ODEs WITH saing all intermediate steps.
    //Shooting_integrator_save_intermediate(ODEdata, N_steps, R_fermi, M_total, inits, mu, lambda);
    //save_data_ode(ODEdata, N_steps, "ODE_data_full.txt");
    // print the results for now (later save them into a txt file).
    //std::cout << "Star with rho_c = " << rho_c << ": radius = " << R_fermi << " [M], mass = " << M_total << " [M_sun]" << std::endl;

    // try the new tabulated EOS system:
    EoStable myDD2("DD2_eos.table");
    sp.EOS = std::make_shared<PolytropicEoS>();

    //vector inits(1.0, 1.0, 1e-20, 1e-20, 100*std::pow(10*rho_c, 2.) );
    //Shooting_integrator_nosave_intermediate(R_fermi, M_total, inits, myDD2, mu, lambda);


    // uncomment this if you want to compute a full MR-diagram

    std::cout << "start loop" << std::endl;

    // create a simple MR-diagram:
    for (unsigned i = 0; i < Nstars; ++i) {

        double rho_start = (i+1)*rho_c;

        // a, alpha, Phi, Psi, P(rho)
        vector inits( {1.0, 1.0, 1e-20, 1e-20, sp.EOS->get_P_from_rho(rho_start)});
        Shooting_integrator_nosave_intermediate(R_fermi, M_total, inits, &sp);
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
