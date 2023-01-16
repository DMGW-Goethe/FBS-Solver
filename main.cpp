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



int Example_Star() {

    double mu = 10.;
    double lambda = 10.;

    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    //auto Polytrope = std::make_shared<PolytropicEoS>();

    double rho_0 = 0.004;
    double phi_0 = 0.06;
    std::vector<integrator::step> steps;

    // define star parameters
    FermionBosonStar fbs(EOS_DD2, mu, lambda, 0., rho_0, phi_0);

    // find omega and integrate
    double omega_0 = 1., omega_1 = 10.;
    //int bisection(double omega_0, double omega_1, int n_mode=0, int max_step=500, double delta_omega=1e-15, int verbose=0);
    fbs.bisection(omega_0, omega_1, 0, 500, 1e-15, 1);
    fbs.evaluate_model(steps, integrator::IntegrationOptions(), "test/fbs.txt");

    std::cout << fbs << std::endl;

    // construct TLN instance from previous instance
    FermionBosonStarTLN fbstln(fbs);

    // find phi_1_0 and integrate
    double phi_1_0_l = phi_0*1e-3, phi_1_0_r = 1e5*phi_0;
    fbstln.bisection_phi_1(phi_1_0_l, phi_1_0_r);
    fbstln.evaluate_model(steps, "test/fbstln.txt");

    std::cout << fbstln << std::endl;


    #ifdef DEBUG_PLOTTING
    //[> see https://github.com/lava/matplotlib-cpp/issues/268 <]
    matplotlibcpp::detail::_interpreter::kill();
    #endif

    return 0;
}


void fillValuesPowerLaw(const double minValue, const double maxValue, std::vector<double>& values, const int power)
{
    if(power == 1)
    {
        const double dValue = double(maxValue - minValue) / double(values.size() - 1);
        for(size_t i = 0; i < values.size(); i++)
            values[i] = minValue + dValue * i;

        return;
    }

    fillValuesPowerLaw(0.0, 1.0, values, 1);

    for(size_t i = 0; i < values.size(); i++)
    {
        values[i] = pow(values[i], power);
        values[i] *= maxValue - minValue;
        values[i] += minValue;
    }
}

void fillValuesLogarithmic(const double minValue, const double maxValue, std::vector<double>& values)
{
    fillValuesPowerLaw(log(minValue), log(maxValue), values, 1);

    for(size_t i = 0; i < values.size(); i++)
        values[i] = exp(values[i]);
}

int main() {
    /*[> see https://github.com/lava/matplotlib-cpp/issues/268
      if this doesn't work, look at the end of the function
    //matplotlibcpp::backend("TkAgg");
    */

    return Example_Star();

    // ----------------------------------------------------------------
    // generate MR curves:
    const unsigned Nstars = 10;     // number of stars in MR curve of constant Phi_c
    const unsigned NstarsPhi = 10;   // number of MR curves of constant rho_c
    const unsigned NstarsNbNf = 2;  // number of MR curves of constand NbNf ratio

    // define some global values:
    double mu = 1.0;        // DM mass
    double lambda = 0.0;    //self-interaction parameter

    constexpr bool calcTln = false;

    // declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    auto Polytrope = std::make_shared<PolytropicEoS>();


    // declare initial conditions:
    double rho_cmin = 1e-8;   // central density of first star (good for DD2 is 0.0005)
    double phi_cmin = 1e-8;    // central value of scalar field of first star
    double rho_cmax = 5e-4;//0.0035;//0.004;
    double phi_cmax = 6e-3;//0.09;//0.07;//0.10;//0.050420813862;

    std::vector<double> rho_c_grid(Nstars, 0.0), phi_c_grid(NstarsPhi, 0.0), NbNf_grid;

    fillValuesPowerLaw(phi_cmin, phi_cmax, phi_c_grid, 1);
    fillValuesPowerLaw(rho_cmin, rho_cmax, rho_c_grid, 1);
    //fillValuesLogarithmic(phi_cmin, phi_cmax, phi_c_grid);
    //fillValuesLogarithmic(rho_cmin, rho_cmax, rho_c_grid);

    //test_EOS(mu, lambda, EOS_DD2, rho_c_grid, phi_c_grid, "plots/DD2_MR_MRphi-plot4.txt");
	// setup to compute a full NS configuration, including tidal deformability:
	std::vector<FermionBosonStar> MRphi_curve;
	std::vector<FermionBosonStarTLN> MRphi_tln_curve;

	// calc the unperturbed equilibrium solutions:
    calc_rhophi_curves(mu, lambda, EOS_DD2, rho_c_grid, phi_c_grid, MRphi_curve);
	// calc the perturbed solutions to get the tidal love number:
    if(calcTln)
	    calc_MRphik2_curve(MRphi_curve, MRphi_tln_curve); // compute the perturbed solutions for TLN
	// save the results in a txt file:

    if(calcTln)
        write_MRphi_curve<FermionBosonStarTLN>(MRphi_tln_curve, "plots/tlncurve_mu1_lambda0_40x40_pow3spacing.txt");
    else
    	write_MRphi_curve<FermionBosonStar>(MRphi_curve, "plots/test.txt");

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
