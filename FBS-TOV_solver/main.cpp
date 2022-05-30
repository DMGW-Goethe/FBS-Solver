#include <iostream> // for i/o e.g. std::cout etc.
#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <fstream>	// file streams
#include <memory>
#include <chrono>   // for timing functionalities
// #include <mpi.h>  // to use MPI (for parallelization) -> needs special compile rules!

#include "vector.hpp"    // include custom 5-vector class
#include "integrator.hpp"
#include "eos.hpp" // include eos container class
#include "nsmodel.hpp"
#include "plotting.hpp"     // to use python/matplotlib inside of c++

// --------------------------------------------------------------------
using clock_type = std::chrono::steady_clock;
using second_type = std::chrono::duration<double, std::ratio<1> >;
typedef std::chrono::time_point<clock_type> time_point;


int main() {
    /* see https://github.com/lava/matplotlib-cpp/issues/268
     * if this doesn't work, look at the end of the function*/
    //matplotlibcpp::backend("TkAgg");

    // for the MR diagram (multiple stars):
    // define some global values:
    // double mu = 1.;        // DM mass
    // double lambda = 0.0;    //self-interaction parameter

    // double rho_c = 0.0002;   // central density
    // double phi_c = 1e-2;    // central value of scalar field
    // double omega_0 = 1., omega_1 =10.;

    // integrate ONE star and save all intermediate values into a txt file:
    // a, alpha, Phi, Psi, P(rho)
    //vector inits(1.0, 1.0, 1e-20, 1e-20, 100*std::pow(10*rho_c, 2.) );
    // solve the ODEs WITH saing all intermediate steps.
    // print the results for now (later save them into a txt file).
    //std::cout << "Star with rho_c = " << rho_c << ": radius = " << R_fermi << " [M], mass = " << M_total << " [M_sun]" << std::endl;

    // try the new tabulated EOS system:
    //auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    //auto EOS_poly = std::make_shared<PolytropicEoS>();

    // declare one FBS object with corresponding initial conditions:
    //FermionBosonStar myFBS(EOS_DD2, mu, lambda, 0.);
    /*
    myFBS.set_initial_conditions(0., rho_c, phi_c);

    // start the bisection search for the correct omega-value in the range [omega_0,omega_1]
    double omega_0 = 1., omega_1 =10.;
    time_point start{clock_type::now()};
    myFBS.bisection(omega_0, omega_1);

    time_point end{clock_type::now()};
    std::cout << "bisection took " << std::chrono::duration_cast<second_type>(end-start).count() << "s" << std::endl;

    time_point start2{clock_type::now()};
    myFBS.evaluate_model("plots/model.txt");
    time_point end2{clock_type::now()};
    std::cout << "evaluation took " << std::chrono::duration_cast<second_type>(end2-start2).count() << "s" << std::endl;
    */

    // ----------------------------------------------------------------
    // generate MR curves:
    const unsigned Nstars = 20;     // number of stars in MR curve of constant Phi_c
    const unsigned NstarsPhi = 1;   // number of MR curves of constant Phi_c
    std::vector<std::vector<double>> myMRcurve(Nstars*NstarsPhi, std::vector<double>(9, 0.0));

    // define some global values:
    double mu = 1.;        // DM mass
    double lambda = 0.0;    //self-interaction parameter

    double rho_c = 0.0002;   // central density of first star
    double phi_c = 1e-20;    // central value of scalar field of first star
    double omega_0 = 1., omega_1 = 10.;  // upper and lower bound for omega in the bisection search

    // declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    auto EOS_poly = std::make_shared<PolytropicEoS>();
    // declare one FBS object with corresponding initial conditions and EOS:
    FermionBosonStar myFBS(EOS_DD2, mu, lambda, 0.);

    time_point start3{clock_type::now()};
    for (unsigned i = 0; i < Nstars; ++i) {
        for (unsigned j = 0; j < NstarsPhi; ++j) {

            double rho_start = i*1e-4 + rho_c;
            double Phi_start = j*0.005 + phi_c;

            // set init data for each star:
            myFBS.set_initial_conditions(0., rho_start, Phi_start);
            myFBS.bisection(omega_0, omega_1);  // compute bisection
            myFBS.evaluate_model("");   // evaluate the model but do not save the intermediate data into txt file

            // extract the relevant values from the solution and fill the MR curve array:
            unsigned index = i + j * Nstars;
            myMRcurve[index][0] = myFBS.M_T;                // total mass
            myMRcurve[index][1] = rho_start;                // central density
            myMRcurve[index][2] = Phi_start;                // central scalar field
            myMRcurve[index][3] = myFBS.R_F*1.476625061;    // fermionic radius
            myMRcurve[index][4] = myFBS.R_F_0*1.476625061;  // fermionic radius where P(r)=0
            myMRcurve[index][5] = myFBS.N_F;                // number of fermions
            myMRcurve[index][6] = myFBS.R_B*1.476625061;    // bosonic radius
            myMRcurve[index][7] = myFBS.N_B;                // number of bosons
            myMRcurve[index][8] = myFBS.N_B / myFBS.N_F;    // ratio N_B / N_F
        }
    }
    time_point end3{clock_type::now()};
    std::cout << "evaluation of "<< Nstars*NstarsPhi <<" stars took " << std::chrono::duration_cast<second_type>(end3-start3).count() << "s" << std::endl;
    std::cout << "average time per evaluation: " << (std::chrono::duration_cast<second_type>(end3-start3).count()/(Nstars*NstarsPhi)) << "s" << std::endl;

    plotting::save_MR_data(myMRcurve, {0,1,2,3,4,5,6,7,9}, {"M","rho_c","phi_c","R_F","R_F_0","N_F","R_B","N_B","N_B/N_F"}, "plots/DD2_MR_MRphi-plot1.txt");
    // ----------------------------------------------------------------

    #ifdef DEBUG_PLOTTING
    /* see https://github.com/lava/matplotlib-cpp/issues/268 */
    matplotlibcpp::detail::_interpreter::kill();
    #endif
    return 0;
}
