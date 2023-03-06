#include <iostream> // for i/o e.g. std::cout etc.
#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
//#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <memory>   // shared_ptr

#include "vector.hpp"    // include custom 5-vector class
#include "integrator.hpp"
#include "eos.hpp" // include eos container class
#include "nsmodel.hpp"
#include "mr_curves.hpp"
#include "plotting.hpp"     // to use python/matplotlib inside of c++

// --------------------------------------------------------------------



int Example_Star() {

    // define bosonic field parameters
    double mu = 10.;
    double lambda = 0.;

    // define EoS
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    //auto Polytrope = std::make_shared<PolytropicEoS>();

    // define central star parameters
    double rho_0 = 0.004;
    double phi_0 = 0.06;
    std::vector<integrator::step> steps;

    // define FBS instance
    FermionBosonStar fbs(EOS_DD2, mu, lambda, 0., rho_0, phi_0);

    // find omega through bisection
    double omega_0 = 1., omega_1 = 10.;
    //int bisection(double omega_0, double omega_1, int n_mode=0, int max_step=500, double delta_omega=1e-15, int verbose=0);
    if(fbs.bisection(omega_0, omega_1) == -1)
        return 0.;

    // evaluate model
    fbs.evaluate_model(steps);

    std::cout << fbs << std::endl;

    // construct TLN instance from previous instance
    FermionBosonStarTLN fbstln(fbs);

    // find phi_1_0 through bisection and evaluate
    double phi_1_0_l = phi_0*1e-3, phi_1_0_r = 1e5*phi_0;
    time_point p1{clock_type::now()};
    fbstln.bisection_phi_1(phi_1_0_l, phi_1_0_r);

    time_point p2{clock_type::now()};
    fbstln.evaluate_model(steps);

    time_point p3{clock_type::now()};

    std::cout << fbstln << std::endl;
    std::cout << "TLN: bisection " << std::chrono::duration_cast<second_type>(p2-p1).count() << "s, evaluation " << std::chrono::duration_cast<second_type>(p3-p2).count() << "s" << std::endl;


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


int create_MR_curve() {

    const unsigned Nstars = 30;     // number of stars in MR curve of constant Phi_c
    const unsigned NstarsPhi = 30;   // number of MR curves of constant rho_c

    // define common star parameters:
    double mu = 1.;        // DM mass
    double lambda = 0.0;   // self-interaction parameter

    constexpr bool calcTln = true;

    // declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    auto Polytrope = std::make_shared<PolytropicEoS>();


    // declare initial conditions:
    double rho_cmin = 0.;   // central density of first star (good for DD2 is 0.0005)
    double phi_cmin = 0.;    // central value of scalar field of first star
    double rho_cmax = 5e-3;
    double phi_cmax = 0.1;

    std::vector<double> rho_c_grid(Nstars, 0.), phi_c_grid(NstarsPhi, 0.), NbNf_grid;

    fillValuesPowerLaw(phi_cmin, phi_cmax, phi_c_grid, 3);
    fillValuesPowerLaw(rho_cmin, rho_cmax, rho_c_grid, 3);
    //fillValuesLogarithmic(phi_cmin, phi_cmax, phi_c_grid);
    //fillValuesLogarithmic(rho_cmin, rho_cmax, rho_c_grid);

	// setup to compute a full NS configuration, including tidal deformability:
	std::vector<FermionBosonStar> MRphi_curve;
	std::vector<FermionBosonStarTLN> MRphi_tln_curve;

	// calc the unperturbed equilibrium solutions
    // this benefits greatly from multiple threads
    calc_rhophi_curves(mu, lambda, EOS_DD2, rho_c_grid, phi_c_grid, MRphi_curve, 2);
	// calc the perturbed solutions to get the tidal love number
    if(calcTln)
	    calc_MRphik2_curve(MRphi_curve, MRphi_tln_curve); // compute the perturbed solutions for TLN
    if(calcTln)
        write_MRphi_curve<FermionBosonStarTLN>(MRphi_tln_curve, "plots/tlncurve_mu1.0_lambda0.0_30x30.e.txt");
    else
        write_MRphi_curve<FermionBosonStar>(MRphi_curve, "plots/curve_mu1.0_lambda0.0_30x30.txt");

    return 0.;
}

int main() {

    // integrate a single star
    // Example_Star();

    // create an MR curve
    create_MR_curve();

    #ifdef DEBUG_PLOTTING
    //[> see https://github.com/lava/matplotlib-cpp/issues/268 <]
    matplotlibcpp::detail::_interpreter::kill();
    #endif
    return 0;
}
