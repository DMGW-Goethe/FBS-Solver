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


void calculate_MRphi_curve(const FermionBosonStar& fbs_model, const std::vector<double>& rho_c_grid, const std::vector<double>& phi_c_grid, std::vector<FermionBosonStar>& MRphi_curve) {

    double omega_0 = 1., omega_1 = 10.;  // upper and lower bound for omega in the bisection search
    MRphi_curve.clear();
    MRphi_curve.reserve(rho_c_grid.size()*phi_c_grid.size());

    time_point start3{clock_type::now()};
    for(unsigned int i = 0; i < rho_c_grid.size(); i++) {
        for(unsigned int j = 0; j < phi_c_grid.size(); j++) {
            FermionBosonStar fbs(fbs_model);
            fbs.set_initial_conditions(0., rho_c_grid[i], phi_c_grid[j]);
            fbs.bisection(omega_0, omega_1);  // compute bisection
            fbs.evaluate_model();   // evaluate the model but do not save the intermediate data into txt file
            MRphi_curve.push_back(fbs);
        }
    }
    time_point end3{clock_type::now()};
    std::cout << "evaluation of "<< MRphi_curve.size() <<" stars took " << std::chrono::duration_cast<second_type>(end3-start3).count() << "s" << std::endl;
    std::cout << "average time per evaluation: " << (std::chrono::duration_cast<second_type>(end3-start3).count()/(MRphi_curve.size())) << "s" << std::endl;

}

void test_EOS(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& phi_c_grid, std::string filename) {

    FermionBosonStar myFBS(EOS, mu, lambda, 0.);    // create model for star
    std::vector<FermionBosonStar> MRphi_curve;

    calculate_MRphi_curve(myFBS, rho_c_grid, phi_c_grid, MRphi_curve); // calculate curves

    if(filename.empty())    // write to file
        return;
    std::ofstream img;
	img.open(filename);
    std::vector<std::string> labels({"M","rho_c","phi_c","R_F","R_F_0","N_F","R_B","N_B","N_B/N_F"});

	if(img.is_open()) {
        // print the labels in line 1:
        for(int i = 0; i < labels.size(); i++)
            img << labels[i] << "\t";
		img << std::endl;

        // print all the data:
        for(auto it = MRphi_curve.begin(); it != MRphi_curve.end(); ++it)
            img << *it << std::endl;
	}
	img.close();
}

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

    // define some global values:
    double mu = 1.;        // DM mass
    double lambda = 0.0;    //self-interaction parameter


    // declare different EOS types:
    auto EOS_DD2 = std::make_shared<EoStable>("EOS_tables/eos_HS_DD2_with_electrons.beta");

    // declare initial conditions:
    double rho_c = 0.0002;   // central density of first star
    double phi_c = 1e-20;    // central value of scalar field of first star

    std::vector<double> rho_c_grid, phi_c_grid;
    for (unsigned i = 0; i < Nstars; ++i)
            rho_c_grid.push_back(i*1e-4 + rho_c);
    for (unsigned j = 0; j < NstarsPhi; ++j)
            phi_c_grid.push_back(j*0.005 + phi_c);


    test_EOS(mu, lambda, EOS_DD2, rho_c_grid, phi_c_grid, "plots/DD2_MR_MRphi-plot2.txt");
    // space for more EOS

    // ----------------------------------------------------------------

    #ifdef DEBUG_PLOTTING
    /* see https://github.com/lava/matplotlib-cpp/issues/268 */
    matplotlibcpp::detail::_interpreter::kill();
    #endif
    return 0;
}
