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

using clock_type = std::chrono::steady_clock;
using second_type = std::chrono::duration<double, std::ratio<1> >;
typedef std::chrono::time_point<clock_type> time_point;



void write_MRphi_curve(const std::vector<FermionBosonStar>& MRphi_curve, std::string filename);

void calc_rhophi_curves(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& phi_c_grid, std::vector<FermionBosonStar>& MRphi_curve);

void calc_NbNf_curves(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& NbNf_grid, std::vector<FermionBosonStar>& MRphi_curve);


#endif
