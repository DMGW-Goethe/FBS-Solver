#ifndef MR_CURVES_HPP
#define MR_CURVES_HPP

#include <iostream> // for i/o e.g. std::cout etc.
#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <fstream>	// file streams
#include <memory>
#include <chrono>   // for timing functionalities
#include <omp.h>

#include "vector.hpp"    // include custom 5-vector class
#include "integrator.hpp"
#include "eos.hpp" // include eos container class
#include "nsmodel.hpp"
#include "fbs_twofluid.hpp"
#include "fbs.hpp"
#include "fbs_tln.hpp"
#include "fps.hpp"

using clock_type = std::chrono::steady_clock;
using second_type = std::chrono::duration<double, std::ratio<1> >;
typedef std::chrono::time_point<clock_type> time_point;


// function to write the FBS results into a txt file. Template functions must be defined in header file:
template <typename T>
void write_MRphi_curve(const std::vector<T>& MRphi_curve, std::string filename) {

    std::ofstream img;
	img.open(filename);
    std::vector<std::string> labels = MRphi_curve.at(0).labels();

	if(img.is_open()) {
        // print the labels in line 1:
        img << "# "; // hashtag so that python recognizes it as a commented line
        for(unsigned i = 0; i < labels.size(); i++)
            img << labels[i] << "\t ";
		img << std::endl;

        // print all the data:
        for(auto it = MRphi_curve.begin(); it != MRphi_curve.end(); ++it)
            img << std::scientific << std::setprecision(10) << *it << std::endl;
	}
	img.close();
}

// function that calculates all FBS/FPS in a rho-phi grid. Adapted to be able to use both FBS and FPS
template <typename U>
void calc_rhophi_curves(std::vector<U>& MRphi_curve, int verbose = 1, int mode = 0) {

    const double omega_0 = 1., omega_1 = 10.;  // upper and lower bound for omega in the bisection search

    unsigned int done = 0;

    time_point start3{clock_type::now()};
    #pragma omp parallel for schedule(dynamic)
    for(unsigned int i = 0; i < MRphi_curve.size(); i++) {
        int bisection_success = MRphi_curve[i].bisection(omega_0, omega_1);  // compute bisection
		if (bisection_success == -1)
            std::cout << "Bisection failed with omega_0=" << omega_0 << ", omega_1=" << omega_1 << " for " << MRphi_curve[i] << std::endl;
        else
            MRphi_curve[i].evaluate_model();   // evaluate the model but do not save the intermediate data into txt file

        #pragma omp atomic
        done++;

        if(verbose >1)
            std::cout << "Progress: "<< float(done) / MRphi_curve.size() * 100.0 << "%" << std::endl;
    }
    time_point end3{clock_type::now()};
    if(verbose > 0) {
        std::cout << "evaluation of "<< MRphi_curve.size() <<" stars took " << std::chrono::duration_cast<second_type>(end3-start3).count() << "s" << std::endl;
        std::cout << "average time per evaluation: " << (std::chrono::duration_cast<second_type>(end3-start3).count()/(MRphi_curve.size())) << "s" << std::endl;
    }

}

//void calc_rhophi_curves(std::vector<FermionBosonStar>& MRphi_curve, int verbose=1);
void calc_rhophi_curves(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& phi_c_grid, std::vector<FermionBosonStar>& MRphi_curve, int verbose=1);

void calc_NbNf_curves(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& NbNf_grid, std::vector<FermionBosonStar>& MRphi_curve);

void calc_MRphik2_curve(const std::vector<FermionBosonStar>& MRphi_curve,  std::vector<FermionBosonStarTLN>& MRphik2_curve, int verbose=1);

void calc_twofluidFBS_curves(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2, const std::vector<double>& rho1_c_grid, const std::vector<double>& rho2_c_grid, std::vector<TwoFluidFBS>& MRphi_curve, double mu=1, double lambda=1);

void calc_rhophi_curves_FPS(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& E_c_grid, std::vector<FermionProcaStar>& MRphi_curve, int verbose=1, int mode = 0);

#endif
