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
    double mu = 1.;        // DM mass
    double lambda = 0.0;    //self-interaction parameter

    double rho_c = 0.002;   // central density
    double phi_c = 1e-2;    // central value of scalar field

    // integrate ONE star and save all intermediate values into a txt file:
    // a, alpha, Phi, Psi, P(rho)
    //vector inits(1.0, 1.0, 1e-20, 1e-20, 100*std::pow(10*rho_c, 2.) );
    // solve the ODEs WITH saing all intermediate steps.
    // print the results for now (later save them into a txt file).
    //std::cout << "Star with rho_c = " << rho_c << ": radius = " << R_fermi << " [M], mass = " << M_total << " [M_sun]" << std::endl;

    // try the new tabulated EOS system:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    auto EOS_poly = std::make_shared<PolytropicEoS>();

    // declare one FBS object with corresponding initial conditions:
    FermionBosonStar myFBS(EOS_DD2, mu, lambda, 0.);
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

    #ifdef DEBUG_PLOTTING
    /* see https://github.com/lava/matplotlib-cpp/issues/268 */
    matplotlibcpp::detail::_interpreter::kill();
    #endif
    return 0;
}
