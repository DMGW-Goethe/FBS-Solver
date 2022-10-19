#include <iostream> // for i/o e.g. std::cout etc.
#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <fstream>	// file streams
#include <memory>

#include "vector.hpp"    // include custom 5-vector class
#include "integrator.hpp"
#include "eos.hpp" // include eos container class
#include "nsmodel.hpp"
#include "mr_curves.hpp"
#include "plotting.hpp"     // to use python/matplotlib inside of c++

// --------------------------------------------------------------------



int test_FBSTLN() {

    double mu = 1.;
    double lambda = 0.;

    //auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    auto Polytrope = std::make_shared<PolytropicEoS>();

    double rho_0 = 1e-20;
    double phi_0 = 0.02;
    std::vector<integrator::step> steps;

    FermionBosonStar fbs(Polytrope, mu, lambda, 0.);
    fbs.set_initial_conditions(rho_0, phi_0);

    double omega_0 = 1., omega_1 = 10.;
    fbs.bisection(omega_0, omega_1);
    fbs.evaluate_model(steps, "test/fbs.txt");

    double H_0 = 1.;
    FermionBosonStarTLN fbstln(fbs);
    fbstln.set_initial_conditions(0., H_0);

    double phi_1_0 = 1e-20, phi_1_1 = 1e3;
    fbstln.bisection_phi_1(phi_1_0, phi_1_1);

    fbstln.evaluate_model(steps, "test/fbstln.txt");

    std::cout << fbs << "\n"
              << fbstln << std::endl;


    #ifdef DEBUG_PLOTTING
    //[> see https://github.com/lava/matplotlib-cpp/issues/268 <]
    matplotlibcpp::detail::_interpreter::kill();
    #endif

    return 0;
}

int main() {
    /*[> see https://github.com/lava/matplotlib-cpp/issues/268
      if this doesn't work, look at the end of the function
    //matplotlibcpp::backend("TkAgg");
    */

    //return test_FBSTLN();

    // ----------------------------------------------------------------
    // generate MR curves:
    const unsigned Nstars = 40;     // number of stars in MR curve of constant Phi_c
    const unsigned NstarsPhi = 40;   // number of MR curves of constant rho_c
    const unsigned NstarsNbNf = 2;  // number of MR curves of constand NbNf ratio

    // define some global values:
    double mu = 2.0;        // DM mass
    double lambda = 0.0;    //self-interaction parameter


    // declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    auto Polytrope = std::make_shared<PolytropicEoS>();

    // declare initial conditions:
    double rho_cmin = 0.0001;   // central density of first star (good for DD2 is 0.0005)
    double phi_cmin = 1e-10;    // central value of scalar field of first star
    double rho_cmax = 0.004;
    double phi_cmax = 0.10;

    double drho = (rho_cmax - rho_cmin) / (Nstars -1.);
    double dphi = (phi_cmax - phi_cmin) / (NstarsPhi -1.);

    std::vector<double> rho_c_grid, phi_c_grid, NbNf_grid;
    for (unsigned i = 0; i < Nstars; ++i) {
            rho_c_grid.push_back(i*drho + rho_cmin);
            std::cout << rho_c_grid[i] << std::endl;}
    for (unsigned j = 0; j < NstarsPhi; ++j) {
            phi_c_grid.push_back(j*dphi + phi_cmin);
            std::cout << phi_c_grid[j] << std::endl;}
    for (unsigned k = 0; k < NstarsNbNf; ++k) {
            NbNf_grid.push_back(k*0.1 + 0.1);
            std::cout << NbNf_grid[k] << std::endl; }

    //test_EOS(mu, lambda, EOS_DD2, rho_c_grid, phi_c_grid, "plots/DD2_MR_MRphi-plot4.txt");
    
	// setup to compute a full NS configuration, including tidal deformability:
	std::vector<FermionBosonStar> MRphi_curve;
	std::vector<FermionBosonStarTLN> MRphi_tln_curve;

	// calc the unperturbed equilibrium solutions:
    calc_rhophi_curves(mu, lambda, EOS_DD2, rho_c_grid, phi_c_grid, MRphi_curve);
	// calc the perturbed solutions to get the tidal love number:
	calc_MRphik2_curve(MRphi_curve, MRphi_tln_curve); // compute the perturbed solutions for TLN
	// save the results in a txt file:
	write_MRphi_curve<FermionBosonStarTLN>(MRphi_tln_curve, "plots/tlncurve_test3dd2.txt");

    // space for more EOS

    // method for the bisection with respect to Nb/Nf:
    //double omega_0 = 1., omega_1 = 10.;
    //FermionBosonStar myFBS(EOS_DD2, mu, lambda, 0.);
    //myFBS.set_initial_conditions(rho_c, phi_c);
    //myFBS.shooting_NbNf_ratio(0.2, 1e-3, omega_0, omega_1); // evaluate model is included
    // myFBS.evaluate_model();

    // calc three MR-curves with different Nb/Nf Ratios!
    //calc_NbNf_curves(mu, lambda, EOS_DD2, rho_c_grid, NbNf_grid, "plots/NbNf_test1.txt");


    // ----------------------------------------------------------------

    #ifdef DEBUG_PLOTTING
    //[> see https://github.com/lava/matplotlib-cpp/issues/268 <]
    matplotlibcpp::detail::_interpreter::kill();
    #endif
    return 0;
}
