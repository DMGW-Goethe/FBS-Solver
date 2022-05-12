#include <iostream> // for i/o e.g. std::cout etc.
#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <fstream>	// file streams
#include <memory>
#include <chrono>   // for timing functionalities
// #include <mpi.h>  // to use MPI (for parallelization) -> needs special compile rules!

#include "vector.hpp"    // include custom 5-vector class
#include "RK45.hpp"
#include "eos.hpp" // include eos container class
#include "nsmodel.hpp"
#include "plotting.hpp"     // to use python/matplotlib inside of c++

// --------------------------------------------------------------------
using clock_type = std::chrono::steady_clock;
using second_type = std::chrono::duration<double, std::ratio<1> >;
typedef std::chrono::time_point<clock_type> time_point;


// saves the integation data into an txt file
void save_integration_data(const std::vector<integrator::step>& res, std::string filename) {

	std::ofstream img;
	img.open(filename);

	if(img.is_open()) {
		img << "# r\t     a\t    alpha\t    Phi\t    Psi\t    P" << std::endl;

        for(auto it = res.begin(); it != res.end(); ++it) {
			img << std::fixed << std::setprecision(10) << it->first;    // radius
            for(int i = 0; i < it->second.size(); i++) img << " " << it->second[i]; // the other variables
            img << std::endl;

		}
	}
	img.close();
}


// find the correct omega-value for a given FBS using bisection in the range [omega_0,omega_1]
// args: FermionBosonStar, vector, min omega, max omega
void Bisection(FermionBosonStar* myFBS, const vector& init_vars, double omega_0, double omega_1) {

    // values/parameters for bisection
    double omega_mid;
    int n_roots_0, n_roots_1, n_roots_mid;   // number of roots in Phi(r) (number of roots corresponds to the modes of the scalar field)
    int n_mode = 0;         // number of the mode we want to compute. (mode n has n roots in the scalar field Phi(r))
    double delta_omega = 1e-16;     // accuracy for omega
    int phi_index = 2;      // index of the Phi-component in the solution array
    int n_max_steps = 1000; // max number of steps in the bisection search after which the bisection search will end
    int i = 0;

    // variables regarding the integration
    integrator::IntegrationOptions intOpts;
    double r_init = 1e-10, r_end= 1000;     // initial integration radius and max radius for integration
    // define events to check for during integration and put them inside of a std::vector:
    integrator::Event phi_neg([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[2] < 0.; });
    integrator::Event phi_pos([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[2] > 0.; });
    std::vector<integrator::Event> events = {phi_neg, phi_pos};     // put the events into the event array
    // declare containers to hold the solution of the integration for the upper- (1), lower- (0) and middle (mid) omega
    std::vector<integrator::step> results_0, results_1, results_mid;

    // find initial values for omega min and omega max
    assert(omega_0 < omega_1);  // if the lower omega (omega_0) is larger than the upper omega (omega_1), we stop the code execution

    // set the lower omega and integrate the ODEs:
    myFBS->omega = omega_0;
    int res =  integrator::RKF45(&(myFBS->dy_dt_static), r_init, init_vars, r_end, (void*) myFBS,  results_0,  events, intOpts);
    n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    // set the upper omega and integrate the ODEs:
    myFBS->omega = omega_1;
    res =  integrator::RKF45(&(myFBS->dy_dt_static), r_init, init_vars, r_end, (void*) myFBS,  results_1,  events, intOpts);
    n_roots_1 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    //assert(omega_0 < omega_1);      // redundant (see a few lines above)
    assert(n_roots_0 != n_roots_1);
    assert(n_roots_0 <= n_mode && n_mode <= n_roots_1);
    std::cout << "start with omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;

    // find right number of zero crossings (roots) cossesponding to the number of modes (n-th mode => n roots)
    // iterate until the upper and lower omega produce results with one root difference
    while(n_roots_1 - n_roots_0 > 1) {
        omega_mid = (omega_0 + omega_1)/2.;
        std::cout << "omega_mid = " << omega_mid << " ->";
        myFBS->omega = omega_mid;
        res =  integrator::RKF45(&(myFBS->dy_dt_static), r_init, init_vars, r_end, (void*) myFBS,  results_mid,  events, intOpts);
        n_roots_mid = events[0].steps.size() + events[1].steps.size() -1;   // number of roots is number of - to + crossings plus + to - crossings
        std::cout << " with n_roots = " << n_roots_mid << std::endl;

        if(n_roots_mid == n_roots_0 || n_roots_mid <= n_mode) {
            n_roots_0 = n_roots_mid;
            omega_0 = omega_mid;
            continue;
        }
        if(n_roots_mid == n_roots_1 || n_roots_mid >= n_mode) {
            n_roots_1 = n_roots_mid;
            omega_1 = omega_mid;
            continue;
        }
    }
    std::cout << "found omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;

    // find right behavior at infty ( Phi(r->infty) = 0 )
    int n_inft_0, n_inft_1, n_inft_mid; // store the sign of Phi at infinity (or at the last r-value)
    myFBS->omega = omega_0; intOpts.save_intermediate=true; // save intermediate is only true because we plot it later with python
    res =  integrator::RKF45(&(myFBS->dy_dt_static), r_init, init_vars, r_end, (void*) myFBS,  results_0,  events, intOpts);
    n_inft_0 = results_0[results_0.size()-1].second[phi_index] > 0.;    // save if sign(Phi(inf)) is positive or negative

    myFBS->omega = omega_1;
    res =  integrator::RKF45(&(myFBS->dy_dt_static), r_init, init_vars, r_end, (void*) myFBS,  results_1,  events, intOpts);
    n_inft_1 = results_1[results_1.size()-1].second[phi_index] > 0.;    // save if sign(Phi(inf)) is positive or negative
    std::cout << "start with omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << " and omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;

    // plot the evolution with python (comment this out later)
    plotting::plot_evolution(results_0, events, {2}, { "Phi_0" });
    plotting::plot_evolution(results_1, events, {2}, { "Phi_1" });
    matplotlibcpp::legend();
    matplotlibcpp::save("bisection_int.png"); matplotlibcpp::close();

    intOpts.save_intermediate=false;
    while(omega_1 - omega_0 > delta_omega && i < n_max_steps) { // iterate until accuracy in omega was reached or max number of steps exceeded
        omega_mid = (omega_0 + omega_1)/2.;
        std::cout << "omega_mid = " << omega_mid << " ->";
        myFBS->omega = omega_mid;
        res =  integrator::RKF45(&(myFBS->dy_dt_static), r_init, init_vars, r_end, (void*) myFBS,  results_mid,  events, intOpts);
        n_inft_mid = results_mid[results_mid.size()-1].second[phi_index] > 0.;  // save if sign(Phi(inf)) is positive or negative
        std::cout << " with n_inft= " << n_inft_mid << std::endl;

        i++;
        // compare the signs of Phi at infinity of the omega-upper, -middle and -lower solution
        // when middle and lower sign are equal, we can move omega_0 to omega_mid
        if(n_inft_mid == n_inft_0) {
            n_inft_0 = n_inft_mid;
            omega_0 = omega_mid;
            continue;
        }
        // when middle and upper sign are equal, we can move omega_1 to omega_mid
        if(n_inft_mid == n_inft_1) {
            n_inft_1 = n_inft_mid;
            omega_1 = omega_mid;
            continue;
        }
    }

    std::cout << "found omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << " and omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;
    myFBS->omega = omega_0;

    // integrate again with saving the intermediate steps to plot them (we remove this later)
    myFBS->omega = omega_0; intOpts.save_intermediate=true;
    res =  integrator::RKF45(&(myFBS->dy_dt_static), r_init, init_vars, r_end, (void*) myFBS,  results_0,  events, intOpts);

    myFBS->omega = omega_1;
    res =  integrator::RKF45(&(myFBS->dy_dt_static), r_init, init_vars, r_end, (void*) myFBS,  results_1,  events, intOpts);
    plotting::plot_evolution(results_0, events, {2}, {"Phi_0"});
    plotting::plot_evolution(results_1, events, {2}, {"Phi_1"});
    matplotlibcpp::legend();
    matplotlibcpp::save("bisection_fin.png"); matplotlibcpp::close();
}


void cumtrapz(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& res) {
    assert(x.size() == y.size() && x.size() > 0);
    res = std::vector<double>(x.size(), 0.);

    for(int i = 1; i < x.size(); i++) {
        res[i] = (x[i]-x[i-1]) * (y[i] + y[i-1])/2.;
        res[i] += res[i-1];
    }
}

void evaluateModel(FermionBosonStar* myFBS, const vector& init_vars) {

    integrator::IntegrationOptions intOpts;
    intOpts.save_intermediate = true;
    double r_init = 1e-10, r_end= 1000;

    integrator::Event M_converged([](const double r, const double dr, const vector& y, const vector& dy, const void *params) {
                                                                                                        double dM_dr = ((1. - 1./y[0]/y[0])/2. + r*dy[0]/y[0]/y[0]/y[0]);                                                                                                         return  dM_dr < 1e-18 ; },
                                                                                                 true);
    std::vector<integrator::Event> events = {M_converged};
    std::vector<integrator::step> results;

    int res =  integrator::RKF45(&(myFBS->dy_dt_static), r_init, init_vars, r_end, (void*) myFBS,  results,  events, intOpts);
    plotting::plot_evolution(results, events, {0,1,2,3,4}, {"a", "alpha", "Phi", "Psi", "P"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log");
    matplotlibcpp::save("evaluation.png"); matplotlibcpp::close();

    auto last_step = results[results.size()-1];
    double M_T = last_step.first / 2. * (1. - 1./last_step.second[0]/last_step.second[0]);

    // Extract the results and put them into a usable form to calculate N_B, N)F
    std::vector<double> r(results.size()), N_B_integrand(results.size()), N_F_integrand(results.size());
    vector v;
    double rho, eps;

    for(int i = 0; i < results.size(); i++) {
        r[i] = results[i].first;
        v = results[i].second;
        N_B_integrand[i] = v[0] * myFBS->omega *  v[2] * v[2] * r[i] * r[i] / v[1];
        myFBS->EOS->callEOS(rho, eps, std::max(0., v[4]));
        N_F_integrand[i] = v[0] * rho * r[i] * r[i] ;
    }

    // Integrate
    std::vector<double> N_F_integrated, N_B_integrated;
    cumtrapz(r, N_F_integrand, N_F_integrated);
    cumtrapz(r, N_B_integrand, N_B_integrated);

    // Find where 99% of N_B,N_F are reached to get the radii
    double N_F = N_F_integrated[N_F_integrated.size()-1], N_B = N_B_integrated[N_B_integrated.size()-1];

    int i_B=-1, i_F=-1;
    for(int i = 0; i < r.size(); i++) {
        if(i_B < 0) {
            if(N_B_integrated[i] > 0.99 * N_B)
                i_B = i;
        }
        if(i_F < 0) {
            if(N_F_integrated[i] > 0.99 * N_F)
                i_F = i;
        }
        if(i_B > 0 && i_F > 0)
            break;
    }
    double R_B = r[i_B], R_F = r[i_F];

    std::cout << "M_T = " << M_T << ", N_B = " << N_B << ", R_B = " << R_B << ", N_F = " << N_F << "R_F = " << R_F << ", N_B/N_F = " << N_B / N_F << std::endl;

}


int main() {
    /* see https://github.com/lava/matplotlib-cpp/issues/268
     * if this doesn't work, look at the end of the function*/
    //matplotlibcpp::backend("TkAgg");

    // for the MR diagram (multiple stars):
    // define some global values:
    double mu = 1.;        // DM mass
    double lambda = 0.0;    //self-interaction parameter

    double rho_c = 0.002;   // central density
    double phi_c = 0.02;    // central value of scalar field

    // integrate ONE star and save all intermediate values into a txt file:
    // a, alpha, Phi, Psi, P(rho)
    //vector inits(1.0, 1.0, 1e-20, 1e-20, 100*std::pow(10*rho_c, 2.) );
    // solve the ODEs WITH saing all intermediate steps.
    // print the results for now (later save them into a txt file).
    //std::cout << "Star with rho_c = " << rho_c << ": radius = " << R_fermi << " [M], mass = " << M_total << " [M_sun]" << std::endl;

    // try the new tabulated EOS system:
    //auto EOS_DD2 = std::make_shared<EoStable>("DD2_eos.table");
    auto EOS_poly = std::make_shared<PolytropicEoS>();

    // declare one FBS object with corresponding initial conditions:
    FermionBosonStar myFBS(EOS_poly, mu, lambda, 0.);
    vector inits =  myFBS.initial_conditions(0., rho_c, phi_c);

    // start the bisection search for the correct omega-value in the range [omega_0,omega_1]
    double omega_0 = 1., omega_1 =10.;
    time_point start{clock_type::now()};
    Bisection(&myFBS, inits, omega_0, omega_1);

    time_point end{clock_type::now()};
    std::cout << "bisection took " << std::chrono::duration_cast<second_type>(end-start).count() << "s" << std::endl;

    time_point start2{clock_type::now()};
    evaluateModel(&myFBS, inits);
    time_point end2{clock_type::now()};
    std::cout << "evaluation took " << std::chrono::duration_cast<second_type>(end2-start2).count() << "s" << std::endl;

    /* see https://github.com/lava/matplotlib-cpp/issues/268 */
    matplotlibcpp::detail::_interpreter::kill();
    return 0;
}
