#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <iostream> // for i/o e.g. std::cout etc.
#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <fstream>	// file streams
#include <vector>
#include <string>
#include <cassert>

#include "integrator.hpp"
#ifdef DEBUG_PLOTTING
#define WITHOUT_NUMPY
#include "matplotlibcpp.h"
#endif


namespace plotting {

void save_integration_data(const std::vector<integrator::step>& results, std::vector<int> plot_components, std::vector<std::string> labels, std::string filename);

void plot_evolution(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events, std::vector<int> plot_components, std::vector<std::string> labels, std::string filename="");

//write any 2D array (2D std::vector) into a txt file

void save_MR_data(std::vector<std::vector<double>> &array, std::vector<int> plot_components, std::vector<std::string> labels, std::string filename);

}
#endif
