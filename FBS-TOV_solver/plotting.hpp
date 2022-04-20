#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <vector>
#include <string>
#include <cassert>

#define WITHOUT_NUMPY
#include "matplotlibcpp.h"
#include "RK45.hpp"



namespace plotting {


void plot_evolution(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events, std::vector<int> plot_components, std::vector<std::string> labels, std::string filename);


}
#endif
