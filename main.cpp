#include <iostream> // for i/o e.g. std::cout etc.
#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
//#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <memory>   // shared_ptr

#include "vector.hpp"    // include custom 5-vector class
#include "integrator.hpp"
#include "eos.hpp" // include eos container class
#include "nsmodel.hpp"
#include "fbs_twofluid.hpp"
#include "mr_curves.hpp"
#include "plotting.hpp"     // to use python/matplotlib inside of c++
#include "utilities.hpp"
#include "fps.hpp"

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


// compute here e.g. a few configurations with the full system and with the effective EOS for different lambda, to check if they produce the same results.
void test_effectiveEOS_pure_boson_star() {

	double mu = 0.5;
	double Lambda_int = 10.;
	double lambda = Lambda_int*8*M_PI*mu*mu;

	// create the phi_c-rho_c-grid:
	const unsigned NstarsPhi = 10;
	const unsigned NstarsRho = 1;
	std::vector<double> rho_c_grid(NstarsRho, 0.0), phi_c_grid(NstarsPhi, 0.0);

	double rho_cmin = 1e-10;	// = [saturation_density] / 0.15 * 2.886376934e-6 * 939.565379 -> includig conversion factors from saturation density to code units
	double rho_cmax = 5e-3;
	double phi_cmin = 1e-5;
	double phi_cmax = 0.055;
	
	utilities::fillValuesPowerLaw(phi_cmin, phi_cmax, phi_c_grid, 1);
    utilities::fillValuesPowerLaw(rho_cmin, rho_cmax, rho_c_grid, 1);

	// declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
	/*auto EOS_APR = std::make_shared<EoStable>("EOS_tables/eos_SRO_APR_SNA_version.beta");
	auto EOS_KDE0v1 = std::make_shared<EoStable>("EOS_tables/eos_SRO_KDE0v1_SNA_version.beta");
	auto EOS_LNS = std::make_shared<EoStable>("EOS_tables/eos_SRO_LNS_SNA_version.beta");
	auto EOS_FSG = std::make_shared<EoStable>("EOS_tables/eos_HS_FSG_with_electrons.beta");*/
	//auto EOS_DD2 = std::make_shared<PolytropicEoS>();
	auto myEffectiveEOS = std::make_shared<EffectiveBosonicEoS>(mu, lambda);


	// compute full self-consistent system:
	// setup to compute a full NS configuration, including tidal deformability:
	std::vector<FermionBosonStar> MRphi_curve; std::vector<FermionBosonStarTLN> MRphi_tln_curve;
	// calc the unperturbed equilibrium solutions:
    calc_rhophi_curves(mu, lambda, EOS_DD2, rho_c_grid, phi_c_grid, MRphi_curve, 2);
	// calc the perturbed solutions to get the tidal love number:
	calc_MRphik2_curve(MRphi_curve, MRphi_tln_curve, 2); // compute the perturbed solutions for TLN
	// save the results in a txt file:
	std::string plotname = "pureBS/1paperplot-TLN-line_pureBS_fullsys-mu_" + std::to_string(mu) + "_Lambdaint_" + std::to_string(Lambda_int);
	//plotname = "tidal_pureEOS_KDE0v1_fullsys";
	//write_MRphi_curve<FermionBosonStar>(MRphi_curve, "plots/" + plotname + ".txt");
	write_MRphi_curve<FermionBosonStarTLN>(MRphi_tln_curve, "plots/" + plotname + ".txt");

	// compute effective two-fluid model:
	std::vector<TwoFluidFBS> twofluid_MRphi_curve;
	calc_twofluidFBS_curves(EOS_DD2, myEffectiveEOS, rho_c_grid, phi_c_grid, twofluid_MRphi_curve, mu, lambda);	// use the effective EOS
	plotname = "pureBS/1paperplot-TLN-line_pureBS_effsys-mu_" + std::to_string(mu) + "_Lambdaint_" + std::to_string(Lambda_int);
	//plotname = "tidal_pureEOS_APR_twofluid";
	//plotname = "rhophi-diagram_DD2_twofluid_" + std::to_string(mu) + "_" + std::to_string(lambda);
	write_MRphi_curve<TwoFluidFBS>(twofluid_MRphi_curve, "plots/" + plotname + ".txt");
}


void test_single_fermion_proca_star() {

	double mu = 1.0;
	double Lambda_int = 0.;
	double lambda = Lambda_int*8*M_PI*mu*mu; //Lambda_int*mu*mu; when using units of Minamisuji (2018)

	double rho_0 = 0.0; //5e-3;
	double E_0 = 0.1 / 1.41421356237309 ;/// 8. / M_PI;

	auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
	FermionProcaStar fps_model(EOS_DD2, mu, lambda, 0.0, rho_0, E_0);    // create model for star

	const double omega_0 = 1., omega_1 = 20.;  // upper and lower bound for omega in the bisection search

	int bisection_success = fps_model.bisection(omega_0, omega_1,0);  // compute bisection
		if (bisection_success == -1) {
            std::cout << "Bisection failed with omega_0=" << omega_0 << ", omega_1=" << omega_1 << std::endl;
		}
		else {
			std::vector<integrator::step> results;
			integrator::IntegrationOptions intOpts = integrator::IntegrationOptions();
            fps_model.evaluate_model(results, intOpts, "plots/FPS/first_fps_test.txt");   // evaluate the model but do not save the intermediate data into txt file
		}

	std::cout << "calculation complete!" << std::endl;
	std::cout << "global quantities:" << std::endl;
	std::cout << std::fixed << std::setprecision(10) << fps_model << std::endl;
}


void compute_fermionProcaStar_grid() {

	double mu = 1.0;
	double Lambda_int = 1.0 * 2.; // factor 2 to convert to the number of minamisuji
	double lambda = Lambda_int*mu*mu;// Lambda_int*8*M_PI*mu*mu;

	// create the E_c-rho_c-grid:
	const unsigned NstarsVec = 60;	// number of pure FPS
	const unsigned NstarsRho = 1;	// number of pure NS
	std::vector<double> rho_c_grid(NstarsRho, 0.0), E_c_grid(NstarsVec, 0.0);

	double rho_cmin = 0.0;
	double rho_cmax = 1e-10;
	double E_cmin = 1e-6;

	double E_cmax = 3.2 / 1.41421356237309; // beware: E_cmax MUST have: E_cmax < m / sqrt(lambda) = 1/4*sqrt(pi*Lambda_int)
	if (lambda > 0. && E_cmax > mu/std::sqrt(lambda)) {
		E_cmax = 0.99*mu/std::sqrt(lambda);
	}
		
	
	utilities::fillValuesPowerLaw(E_cmin, E_cmax, E_c_grid, 1);
    utilities::fillValuesPowerLaw(rho_cmin, rho_cmax, rho_c_grid, 1);

	// declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");

	std::vector<FermionProcaStar> FPS_curve;

	calc_rhophi_curves_FPS(mu, lambda, EOS_DD2, rho_c_grid, E_c_grid, FPS_curve, 2, 0);

	std::string plotname = "FPS/Minamisuji_fig3_reproduction-mu_" + std::to_string(mu) + "_Lambdaint_" + std::to_string(Lambda_int);
	write_MRphi_curve<FermionProcaStar>(FPS_curve, "plots/" + plotname + ".txt");
	
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

    utilities::fillValuesPowerLaw(phi_cmin, phi_cmax, phi_c_grid, 3);
    utilities::fillValuesPowerLaw(rho_cmin, rho_cmax, rho_c_grid, 3);
    //utilities::fillValuesLogarithmic(phi_cmin, phi_cmax, phi_c_grid);
    //utilities::fillValuesLogarithmic(rho_cmin, rho_cmax, rho_c_grid);

	// setup to compute a full NS configuration, including tidal deformability:
	std::vector<FermionBosonStar> MRphi_curve;
	std::vector<FermionBosonStarTLN> MRphi_tln_curve;

	// calc the unperturbed equilibrium solutions
    // this benefits greatly from multiple threads
    calc_rhophi_curves(mu, lambda, EOS_DD2, rho_c_grid, phi_c_grid, MRphi_curve, 2);
	// calc the perturbed solutions to get the tidal love number
    if(calcTln)
	    calc_MRphik2_curve(MRphi_curve, MRphi_tln_curve, 2); // compute the perturbed solutions for TLN

    // save results in file
    if(calcTln)
        write_MRphi_curve<FermionBosonStarTLN>(MRphi_tln_curve, "plots/tlncurve_mu1.0_lambda0.0_30x30.txt");
    else
        write_MRphi_curve<FermionBosonStar>(MRphi_curve, "plots/curve_mu1.0_lambda0.0_30x30.txt");

    return 0.;
}

int main() {

    // integrate a single star
    // Example_Star();

    // create an MR curve
    //create_MR_curve();
	// ----------------------------------------------------------------

	// test two-fluid EOS with effective bosonic EoS:
	//test_effectiveEOS_pure_boson_star();
	//test_single_fermion_proca_star();
	compute_fermionProcaStar_grid();

    // ----------------------------------------------------------------

    #ifdef DEBUG_PLOTTING
    //[> see https://github.com/lava/matplotlib-cpp/issues/268 <]
    matplotlibcpp::detail::_interpreter::kill();
    #endif
    return 0;
}
