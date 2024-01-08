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
using second_type = std::chrono::duration<NUMERIC, std::ratio<1> >;
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

    const NUMERIC omega_0 = 1._num, omega_1 = 10._num;  // upper and lower bound for omega in the bisection search

    unsigned int done = 0;

    time_point start3{clock_type::now()};
    #pragma omp parallel for schedule(dynamic)
    for(unsigned int i = 0; i < MRphi_curve.size(); i++) {
        int bisection_success = MRphi_curve[i].bisection(omega_0, omega_1, mode);  // compute bisection
		if (bisection_success == -1)
            std::cout << "Bisection failed with omega_0=" << omega_0 << ", omega_1=" << omega_1 << " for " << MRphi_curve[i] << std::endl;
        else
            MRphi_curve[i].evaluate_model();   // evaluate the model but do not save the intermediate data into txt file

        #pragma omp atomic
        done++;

        if(verbose >1)
            std::cout << "Progress: "<< float(done) / MRphi_curve.size() * 100.0_num << "%" << std::endl;
    }
    time_point end3{clock_type::now()};
    if(verbose > 0) {
        std::cout << "evaluation of "<< MRphi_curve.size() <<" stars took " << std::chrono::duration_cast<second_type>(end3-start3).count() << "s" << std::endl;
        std::cout << "average time per evaluation: " << (std::chrono::duration_cast<second_type>(end3-start3).count()/(MRphi_curve.size())) << "s" << std::endl;
    }

}

// compute curves of constant rho_c and Nb/(Nf+Nb)-ratio:
template <typename V>
void calc_NbNf_curves(NUMERIC mu, NUMERIC lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<NUMERIC>& rho_c_grid, const std::vector<NUMERIC>& NbNf_grid, std::vector<V>& MRphi_curve, int verbose) {
    
	MRphi_curve.clear();
    MRphi_curve.reserve(rho_c_grid.size()*NbNf_grid.size());
    // set initial conditions for every star in the list:
    for (unsigned j = 0; j < NbNf_grid.size(); j++) {
        for(unsigned i = 0; i < rho_c_grid.size(); i++) {
            V fbs(EOS, mu, lambda, 0._num, rho_c_grid[i], 0._num);  // create a star
            MRphi_curve.push_back(fbs);
        }
    }
    // compute the MR-diagrams:
    NUMERIC omega_0 = 1._num, omega_1 = 10._num;
    unsigned int done = 0;
	std::cout << "start loop" << std::endl;
    #pragma omp parallel for schedule(dynamic, 4)
    for(unsigned j = 0; j < rho_c_grid.size() ; j++) {
		for(unsigned i = 0; i < NbNf_grid.size(); i++) {
            int index = i*rho_c_grid.size() + j;
            MRphi_curve[index].shooting_NbNf_ratio(NbNf_grid[i], 1e-4_num, omega_0, omega_1, 0, 200, 1e-16_num, verbose);  // compute star with set Nb/(Nf+Nb) ratio
            MRphi_curve[index].evaluate_model();
			done++;
            std::cout << "Progress: "<< float(done) / (NbNf_grid.size() * rho_c_grid.size()) * 100.0 << "%" << std::endl;
        }
    }
    std::cout << "end loop" << std::endl;
}

void calc_rhophi_curves(NUMERIC mu, NUMERIC lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<NUMERIC>& rho_c_grid, const std::vector<NUMERIC>& phi_c_grid, std::vector<FermionBosonStar>& MRphi_curve, int verbose=1);

void calc_MRphik2_curve(const std::vector<FermionBosonStar>& MRphi_curve,  std::vector<FermionBosonStarTLN>& MRphik2_curve, int verbose=1);

void calc_twofluidFBS_curves(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2, const std::vector<NUMERIC>& rho1_c_grid, const std::vector<NUMERIC>& rho2_c_grid, std::vector<TwoFluidFBS>& MRphi_curve, NUMERIC mu=1._num, NUMERIC lambda=1._num);

void calc_rhophi_curves_FPS(NUMERIC mu, NUMERIC lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<NUMERIC>& rho_c_grid, const std::vector<NUMERIC>& E_c_grid, std::vector<FermionProcaStar>& MRphi_curve, int verbose=1, int mode = 0);

#endif
