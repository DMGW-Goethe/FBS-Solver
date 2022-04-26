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
#include "plotting.hpp"

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


void Bisection(NSmodel*m, const vector& init_vars, double omega_0, double omega_1) {
    // values for bisection
    double omega_mid;
    int n_zc_0, n_zc_1, n_zc_mid;
    int n_mode = 0;
    double delta_omega = 1e-16;
    int phi_index = 2;
    int n_max_steps = 1000;
    int i = 0;

    // variables regarding the integration
    integrator::IntegrationOptions intOpts;
    double r_init = 1e-10, r_end= 1000;
    integrator::Event phi_neg([](const double r, const vector y, const void*params) { return y[2] < 0.; });
    integrator::Event phi_pos([](const double r, const vector y, const void*params) { return y[2] > 0.; });
    std::vector<integrator::Event> events = {phi_neg, phi_pos};
    std::vector<integrator::step> results_0, results_1, results_mid;

    // find initial values
    assert(omega_0 < omega_1);
    m->omega = omega_0;
    int res =  integrator::RKF45(&(m->dy_dt_static), r_init, init_vars, r_end, (void*) m,  results_0,  events, intOpts);
    n_zc_0 = events[0].steps.size() + events[1].steps.size() -1;

    m->omega = omega_1;
    res =  integrator::RKF45(&(m->dy_dt_static), r_init, init_vars, r_end, (void*) m,  results_1,  events, intOpts);
    n_zc_1 = events[0].steps.size() + events[1].steps.size() - 1;
    assert(omega_0 < omega_1);
    assert(n_zc_0 != n_zc_1);
    assert(n_zc_0 <= n_mode && n_mode <= n_zc_1);
    std::cout << "start with omega_0 =" << omega_0 << " with n_zc=" << n_zc_0 << " and omega_1=" << omega_1 << " with n_zc=" << n_zc_1 << std::endl;

    // find right number of zero crossings
    while(n_zc_1 - n_zc_0 > 1) {
        omega_mid = (omega_0 + omega_1)/2.;
        std::cout << "omega_mid = " << omega_mid << " ->";
        m->omega = omega_mid;
        res =  integrator::RKF45(&(m->dy_dt_static), r_init, init_vars, r_end, (void*) m,  results_mid,  events, intOpts);
        n_zc_mid = events[0].steps.size() + events[1].steps.size() -1;
        std::cout << " with n_zc = " << n_zc_mid << std::endl;

        if(n_zc_mid == n_zc_0 || n_zc_mid <= n_mode) {
            n_zc_0 = n_zc_mid;
            omega_0 = omega_mid;
            continue;
        }
        if(n_zc_mid == n_zc_1 || n_zc_mid >= n_mode) {
            n_zc_1 = n_zc_mid;
            omega_1 = omega_mid;
            continue;
        }
    }
    std::cout << "found omega_0 =" << omega_0 << " with n_zc=" << n_zc_0 << " and omega_1=" << omega_1 << " with n_zc=" << n_zc_1 << std::endl;

    // find right behavior at infty
    int n_inft_0, n_inft_1, n_inft_mid;
    m->omega = omega_0; intOpts.save_intermediate=true;
    res =  integrator::RKF45(&(m->dy_dt_static), r_init, init_vars, r_end, (void*) m,  results_0,  events, intOpts);
    n_inft_0 = results_0[results_0.size()-1].second[phi_index] > 0.;

    m->omega = omega_1;
    res =  integrator::RKF45(&(m->dy_dt_static), r_init, init_vars, r_end, (void*) m,  results_1,  events, intOpts);
    n_inft_1 = results_1[results_1.size()-1].second[phi_index] > 0.;
    std::cout << "start with omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << " and omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;

    plotting::plot_evolution(results_0, events, {2}, { "Phi_0" });
    plotting::plot_evolution(results_1, events, {2}, { "Phi_1"});
    matplotlibcpp::legend();
    matplotlibcpp::save("bisection_int.png"); matplotlibcpp::close();

    intOpts.save_intermediate=false;
    while(omega_1 - omega_0 > delta_omega && i < n_max_steps) {
        omega_mid = (omega_0 + omega_1)/2.;
        std::cout << "omega_mid = " << omega_mid << " ->";
        m->omega = omega_mid;
        res =  integrator::RKF45(&(m->dy_dt_static), r_init, init_vars, r_end, (void*) m,  results_mid,  events, intOpts);
        n_inft_mid = results_mid[results_mid.size()-1].second[phi_index] > 0.;
        std::cout << " with n_inft= " << n_inft_mid << std::endl;

        i++;
        if(n_inft_mid == n_inft_0) {
            n_inft_0 = n_inft_mid;
            omega_0 = omega_mid;
            continue;
        }
        if(n_inft_mid == n_inft_1) {
            n_inft_1 = n_inft_mid;
            omega_1 = omega_mid;
            continue;
        }
    }

    std::cout << "found omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << " and omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;
    m->omega = omega_0;

    m->omega = omega_0; intOpts.save_intermediate=true;
    res =  integrator::RKF45(&(m->dy_dt_static), r_init, init_vars, r_end, (void*) m,  results_0,  events, intOpts);

    m->omega = omega_1;
    res =  integrator::RKF45(&(m->dy_dt_static), r_init, init_vars, r_end, (void*) m,  results_1,  events, intOpts);
    plotting::plot_evolution(results_0, events, {2}, {"Phi_0"});
    plotting::plot_evolution(results_1, events, {2}, {"Phi_1"});
    matplotlibcpp::legend();
    matplotlibcpp::save("bisection_fin.png"); matplotlibcpp::close();


}



int main() {
    /* see https://github.com/lava/matplotlib-cpp/issues/268
     * if this doesn't work, look at the end of the function*/
    //matplotlibcpp::backend("TkAgg");

    // for the MR diagram (multiple stars):
    // define some global values:
    double mu = 1.;        // DM mass
    double lambda = 0.0;    //self-interaction parameter

    double rho_c = 0.002; // central density
    double phi_c = 0.02;

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

    vector inits =  m.initial_conditions(0., rho_c, phi_c);
    double omega_0 = 1., omega_1 =10.;
    Bisection(&m, inits, omega_0, omega_1);

    /* see https://github.com/lava/matplotlib-cpp/issues/268 */
    matplotlibcpp::detail::_interpreter::kill();
    return 0;
}
