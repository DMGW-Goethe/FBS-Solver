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
    NUMERIC mu = 10._num;
    NUMERIC lambda = 0._num;

    // define EoS
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    //auto Polytrope = std::make_shared<PolytropicEoS>();

    // define central star parameters
    NUMERIC rho_0 = 0.004_num;
    NUMERIC phi_0 = 0.06_num;
    std::vector<integrator::step> steps;

    // define FBS instance
    FermionBosonStar fbs(EOS_DD2, mu, lambda, 0._num, rho_0, phi_0);

    // find omega through bisection
    NUMERIC omega_0 = 1._num, omega_1 = 10._num;
    //int bisection(NUMERIC omega_0, NUMERIC omega_1, int n_mode=0, int max_step=500, NUMERIC delta_omega=1e-15, int verbose=0);
    if(fbs.bisection(omega_0, omega_1) == -1)
        return 0._num;

    // evaluate model
    fbs.evaluate_model(steps);

    std::cout << fbs << std::endl;

    // construct TLN instance from previous instance
    FermionBosonStarTLN fbstln(fbs);

    // find phi_1_0 through bisection and evaluate
    NUMERIC phi_1_0_l = phi_0*1e-3_num, phi_1_0_r = 1e5_num*phi_0;
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

	NUMERIC mu = 0.5_num;
	NUMERIC Lambda_int = 10._num;
	NUMERIC lambda = Lambda_int*8._num*M_PI*mu*mu;

	// create the phi_c-rho_c-grid:
	const unsigned NstarsPhi = 10;
	const unsigned NstarsRho = 1;
	std::vector<NUMERIC> rho_c_grid(NstarsRho, 0._num), phi_c_grid(NstarsPhi, 0._num);

	NUMERIC rho_cmin = 1e-10_num;	// = [saturation_density] / 0.15 * 2.886376934e-6 * 939.565379 -> includig conversion factors from saturation density to code units
	NUMERIC rho_cmax = 5e-3_num;
	NUMERIC phi_cmin = 1e-5_num;
	NUMERIC phi_cmax = 0.055_num;
	
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

	NUMERIC mu = 1.0_num;
	NUMERIC Lambda_int = 0.0_num;
	NUMERIC lambda = Lambda_int*8._num*M_PI*mu*mu; //Lambda_int*mu*mu; when using units of Minamisuji (2018)

	NUMERIC sat_to_code = 0.15_num * 2.886376934e-6_num * 939.565379_num;	// conversion factor from saturation density to code units
	NUMERIC rho_0 = 9.0_num * sat_to_code; //5e-3;
	NUMERIC E_0 = 0.18_num ;/// 8. / M_PI;

	auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
	FermionProcaStar fps_model(EOS_DD2, mu, lambda, 0._num, rho_0, E_0);    // create model for star
	FermionBosonStar fbs_model(EOS_DD2, mu, lambda, 0._num, rho_0, E_0);	// also create a FBS for comparisons sake

	std::string plotname = "FPS/posterplot1-radial-FPS-5rho_sat-m01-l0-E_0-" + std::to_string(E_0);
	plotname = "FPS/bisection-fail-test-radial-profiles";
	const NUMERIC omega_0 = 1.0_num /*2.12055_num*/, omega_1 = 10.0_num/* 3.2499911682_num */;  // upper and lower bound for omega in the bisection search

	int bisection_success1 = fps_model.bisection(omega_0, omega_1,0);  // compute bisection
	//fps_model.omega = 2.1399911682_num;
		if (bisection_success1 == -5) {
            std::cout << "Bisection failed with omega_0=" << omega_0 << ", omega_1=" << omega_1 << std::endl;
		}
		else {
			std::vector<integrator::step> results1;
			integrator::IntegrationOptions intOpts = integrator::IntegrationOptions();
            fps_model.evaluate_model(results1, intOpts, "plots/"+plotname+".txt");   // evaluate the model and save the intermediate data into txt file
		}
	// do the same for the FBS:
	/*
	plotname = "FPS/posterplot1-radial-FBS-5rho_sat-m01-l0-E_0-" + std::to_string(E_0);
	int bisection_success2 = fbs_model.bisection(omega_0, omega_1,0);  // compute bisection
		if (bisection_success2 == -1) {
            std::cout << "Bisection failed with omega_0=" << omega_0 << ", omega_1=" << omega_1 << std::endl;
		}
		else {
			std::vector<integrator::step> results2;
			integrator::IntegrationOptions intOpts = integrator::IntegrationOptions();
            fbs_model.evaluate_model(results2, intOpts, "plots/"+plotname+".txt");   // evaluate the model and save the intermediate data into txt file
		}*/

	std::cout << "calculation complete!" << std::endl;
	std::cout << "global quantities:" << std::endl;
	std::cout << std::fixed << std::setprecision(10) << fps_model << std::endl;
	std::cout << std::fixed << std::setprecision(10) << fbs_model << std::endl;
}


void compute_fermionProcaStar_grid() {

	NUMERIC mu = 2.0_num;
	NUMERIC Lambda_int = 50.0_num; // factor 2 to convert to the number of minamisuji
	NUMERIC lambda = Lambda_int*8._num*M_PI*mu*mu;// Lambda_int*8*M_PI*mu*mu;

	// create the E_c-rho_c-grid:
	const unsigned NstarsVec = 40;	// number of pure FPS
	const unsigned NstarsRho = 40;	// number of pure NS
	std::vector<NUMERIC> rho_c_grid(NstarsRho, 0.0_num), E_c_grid(NstarsVec, 0.0_num);

	NUMERIC sat_to_code = 0.15_num * 2.886376934e-6_num * 939.565379_num;	// conversion factor from saturation density to code units
	NUMERIC rho_cmin = 1e-4_num * sat_to_code;
	NUMERIC rho_cmax = 10.0_num * sat_to_code;
	NUMERIC E_cmin = 1e-5_num;

	NUMERIC E_cmax = 0.1_num; //3.2 / 1.41421356237309; // beware: E_cmax MUST have: E_cmax < m / sqrt(lambda) = 1/4*sqrt(pi*Lambda_int)
	if (lambda > 0._num && E_cmax > mu/std::sqrt(lambda)) {
		E_cmax = 0.99_num*mu/std::sqrt(lambda);
	}
		
	
	utilities::fillValuesPowerLaw(E_cmin, E_cmax, E_c_grid, 1);
    utilities::fillValuesPowerLaw(rho_cmin, rho_cmax, rho_c_grid, 1);

	// declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
	auto EOS_POLYTROPE = std::make_shared<PolytropicEoS>();

	std::vector<FermionProcaStar> FPS_curve;

	calc_rhophi_curves_FPS(mu, lambda, EOS_DD2, rho_c_grid, E_c_grid, FPS_curve, 2, 0);

	std::string plotname = "FPS/Minamisuji_fig3_reproduction-unitchange-new-mu_" + std::to_string(mu) + "_Lambdaint_" + std::to_string(Lambda_int);
	plotname = "FPS/bisection-fail-test7";
	write_MRphi_curve<FermionProcaStar>(FPS_curve, "plots/" + plotname + ".txt");
	std::cout << "Finished FPS!" << std::endl;

	std::vector<FermionBosonStar> MRphi_curve;
	// calc the unperturbed equilibrium solutions:
    //calc_rhophi_curves(mu, lambda, EOS_DD2, rho_c_grid, E_c_grid, MRphi_curve,2);
	plotname = "FPS/posterplot2-MR-FBS-m1-l0-rhoc-9";
	//write_MRphi_curve<FermionBosonStar>(MRphi_curve, "plots/" + plotname + ".txt");
	
}

int create_MR_curve() {

    const unsigned Nstars = 30;     // number of stars in MR curve of constant Phi_c
    const unsigned NstarsPhi = 30;   // number of MR curves of constant rho_c

    // define common star parameters:
    NUMERIC mu = 1._num;        // DM mass
    NUMERIC lambda = 0._num;   // self-interaction parameter

    constexpr bool calcTln = true;

    // declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
    auto Polytrope = std::make_shared<PolytropicEoS>();


    // declare initial conditions:
    NUMERIC rho_cmin = 0._num;   // central density of first star (good for DD2 is 0.0005)
    NUMERIC phi_cmin = 0._num;    // central value of scalar field of first star
    NUMERIC rho_cmax = 5e-3_num;
    NUMERIC phi_cmax = 0.1_num;

    std::vector<NUMERIC> rho_c_grid(Nstars, 0._num), phi_c_grid(NstarsPhi, 0._num), NbNf_grid;

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


// calculates lines of constnt boson number Nb/(Nf+Nb)
void compute_boson_fraction_grid() {

	NUMERIC mu = 1.0_num;
	NUMERIC Lambda_int = 0.0_num; // factor 2 to convert to the number of minamisuji
	NUMERIC lambda = Lambda_int*8._num*M_PI*mu*mu;// Lambda_int*8*M_PI*mu*mu;

	// create the rho_c-grid to boson-fraction grid:
	const unsigned NstarsRho = 30;	// number of stars per NbNf line
	const unsigned NstarsNbNf = 3;	// number of NbNf lines
	std::vector<NUMERIC> rho_c_grid(NstarsRho, 0.0_num), NbNf_grid(NstarsNbNf, 0.0_num);

	NUMERIC sat_to_code = 0.15_num * 2.886376934e-6_num * 939.565379_num;	// conversion factor from saturation density to code units
	NUMERIC rho_cmin = 1e-0_num * sat_to_code;
	NUMERIC rho_cmax = 9.0_num * sat_to_code;
	NUMERIC NbNf_min = 0.0_num;
	NUMERIC NbNf_max = 0.2_num;

	utilities::fillValuesPowerLaw(NbNf_min, NbNf_max, NbNf_grid, 1);
    utilities::fillValuesPowerLaw(rho_cmin, rho_cmax, rho_c_grid, 2);

	// declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");
	auto EOS_POLYTROPE = std::make_shared<PolytropicEoS>();

	std::vector<FermionBosonStar> FBS_curve;
	calc_NbNf_curves(mu, lambda, EOS_DD2, rho_c_grid, NbNf_grid, FBS_curve, 1);
	std::string plotname = "FPS/test-NbNB-ratio-mu_" + std::to_string(mu) + "_Lambdaint_" + std::to_string(Lambda_int);
	//plotname = "FPS/bisection-fail-test7";
	write_MRphi_curve<FermionBosonStar>(FBS_curve, "plots/" + plotname + ".txt");
	std::cout << "Finished FBS!" << std::endl;
	
		/*
	std::vector<FermionProcaStar> FPS_curve;
    calc_NbNf_curves_FPS(mu, lambda, EOS_DD2, rho_c_grid, NbNf_grid, FPS_curve);
	plotname = "FPS/posterplot2-MR-FBS-m1-l0-rhoc-9";
	write_MRphi_curve<FermionProcaStar>(FPS_curve, "plots/" + plotname + ".txt");
		*/
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
	//compute_fermionProcaStar_grid();
	compute_boson_fraction_grid();

    // ----------------------------------------------------------------

    #ifdef DEBUG_PLOTTING
    //[> see https://github.com/lava/matplotlib-cpp/issues/268 <]
    matplotlibcpp::detail::_interpreter::kill();
    #endif
    return 0;
}
