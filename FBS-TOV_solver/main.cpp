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
#include "nsmodel.hpp"

// --------------------------------------------------------------------


void save_integration_data(const std::vector<integrator::step>& res, std::string filename) {

	std::ofstream img;
	img.open(filename);

	if(img.is_open()) {
		img << "# r\t     a\t    alpha\t    Phi\t    Psi\t    P" << std::endl;

        for(auto it = res.begin(); it != res.end(); ++it) {
			img << std::fixed << std::setprecision(10) << it->first;
            for(int i = 0; i < it->second.size(); i++) img << " " << it->second[i];
            img << std::endl;

		}
	}
	img.close();
}

// integrating the system of equations using the shooting method:
// for shoooting method see e.g. https://sdsawtelle.github.io/blog/output/shooting-method-ode-numerical-solution.html
// parameters: star Radius, Mass, initial values, mu, lambda
// returns only the Mass and fermionic radius of ONE star
void Shooting_integrator_nosave_intermediate(double& Rfermionic, double& Mtot, vector& init_vars, NSmodel* m) {

    // define the evolution variables
    double omega_bar = 1.0;             // initial guess for omega and running variable
    double omega_0bestguess; // this will be the best guess for omega
    double delta_omega = 1e-3;           // to track the quality of the derivative for finite-difference
    const double delta_omega_max = 100.;
    const double root_desired_error = 1e-8;    // desired accuracy for the root F(omega)=0

    const double init_r = 1e-10;        // initial radius (cannot be zero due to coordinate singularity)
    const double r_end = 1000.;
    int max_step = 1000000;

    // declare the variables to be able to use them in the whole function scope
    double my_r = init_r;
    vector first_result, second_result;
    double first_shot, second_shot;

    integrator::Event Rf_condition([](const double r, const vector y, const void*params) { return y[4] < 1e-14; });
    integrator::Event phi_condition([](const double r, const vector y, const void*params) { return y[2] < 1e-14; });
    integrator::Event Rf_phi_condition([](const double r, const vector y, const void*params) { return std::max(y[2], y[4]) < 1e-14; }, true);;
    std::vector<integrator::Event> events = {Rf_condition, phi_condition, Rf_phi_condition};
    std::vector<integrator::step> results;

    // integrate one TOV star until a stable configuration has been reached (including shooting method):
    // outer shooting method loop:
    while(true) {

        // reset the initial values before each ODE itegration (only omega should change each iteration):
        Rfermionic = 0.0;
        Mtot = 0.0;

        // inner ODE integration loop first solution:
        m->omega = omega_bar;
        int res =  integrator::RKF45(&(m->dy_dt_static), init_r, init_vars, r_end, (void*) m, max_step,  results,  events);
        first_result = results[results.size()-1].second; my_r = results[results.size()-1].first;

        Rfermionic = events[0].steps.at(0).first;
        std::cout << "stopped with events[0].active=" << events[0].active << ", events[1].active=" << events[1].active << ", events[2].active=" << events[2].active<< std::endl;
        for(auto it = events.begin(); it != events.end(); ++it)
            it->reset();

        // second solution
        m->omega = omega_bar + delta_omega;
        res =  integrator::RKF45(&(m->dy_dt_static), init_r, init_vars, r_end, (void*) m, max_step,  results,  events);
        for(auto it = events.begin(); it != events.end(); ++it)
            it->reset();
        second_result = results[results.size()-1].second;

        first_shot = first_result[3];
        second_shot = second_result[3];

        // take result from the integrator and try to obtain a new value for omega:
        // Newton-Raphson method to step in the right direction for omega
        std::cout << "first result = " << first_result <<  "second result = " << second_result << std::endl;

        if ( std::abs( first_shot ) < root_desired_error) { // stop the loop once a satisfactory level of accuracy has been found
            break;
		}
        if ( first_shot == second_shot) {
            std::cout << "omega = " << omega_bar << "no change detected. increase delta_omega=" << delta_omega << std::endl;
            delta_omega *= 10.;
            if(delta_omega > delta_omega_max)
                break;
            continue;
        }
        omega_0bestguess =  omega_bar - first_shot * delta_omega / ( second_shot - first_shot );
		omega_bar = omega_0bestguess; // let the best guess become the starting point for the next iteration step

        std::cout << "one_omega_iteration done! omega = " << omega_0bestguess << std::endl;
    }

    m->omega = omega_0bestguess;
    //save_integration_data(results, "test");

    Mtot = 0.5 * my_r * ( 1.0 - 1.0/(std::pow(first_result[0],2) ));
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

    // for the MR diagram (multiple stars):
    const unsigned Nstars = 10;
    double MRarray[Nstars][3] = {0.};

    // define some global values:
    double mu = 0.0;        // DM mass
    double lambda = 0.0;    //self-interaction parameter
    double R_fermi = 0.0; double M_total = 0.0; // fermionic radius of star and total mass of star

    double rho_c = 2e-4; // central density
    double phi_c = 0.;

    // integrate ONE star and save all intermediate values into a txt file:
    // a, alpha, Phi, Psi, P(rho)
    //vector inits(1.0, 1.0, 1e-20, 1e-20, 100*std::pow(10*rho_c, 2.) );
    // solve the ODEs WITH saing all intermediate steps.
    // print the results for now (later save them into a txt file).
    //std::cout << "Star with rho_c = " << rho_c << ": radius = " << R_fermi << " [M], mass = " << M_total << " [M_sun]" << std::endl;

    // try the new tabulated EOS system:
    auto EOS_DD2 = std::make_shared<EoStable>("DD2_eos.table");
    auto EOS_poly = std::make_shared<PolytropicEoS>();

    NSmodel m( EOS_poly, mu, lambda, 0.);

    // uncomment this if you want to compute a full MR-diagram

    std::cout << "start loop" << std::endl;

    // create a simple MR-diagram:
    for (unsigned i = 0; i < Nstars; ++i) {

        double rho_start = (i+1)*rho_c;

        // a, alpha, Phi, Psi, P(rho)
        vector inits =  m.initial_conditions(0., rho_start, phi_c);
        Shooting_integrator_nosave_intermediate(R_fermi, M_total, inits, &m);
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
    //save_data_MR(MRarray, Nstars, "MRtest_DD2eos.txt");


    // print the results for now (later save them into a txt file).
    //std::cout << "Star with rho_c = " << rho_c << ": radius = " << R_fermi << " [M], mass = " << M_total << " [M_sun?]" << std::endl;


    return 0;
}
