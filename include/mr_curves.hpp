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

void calc_rhophi_curves(std::vector<FermionBosonStar>& MRphi_curve, int verbose=1);
void calc_rhophi_curves(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& phi_c_grid, std::vector<FermionBosonStar>& MRphi_curve, int verbose=1);

void calc_NbNf_curves(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& NbNf_grid, std::vector<FermionBosonStar>& MRphi_curve);

void calc_MRphik2_curve(const std::vector<FermionBosonStar>& MRphi_curve,  std::vector<FermionBosonStarTLN>& MRphik2_curve, int verbose=1);

void calc_twofluidFBS_curves(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2, const std::vector<double>& rho1_c_grid, const std::vector<double>& rho2_c_grid, std::vector<TwoFluidFBS>& MRphi_curve, double mu=1, double lambda=1);

#endif
