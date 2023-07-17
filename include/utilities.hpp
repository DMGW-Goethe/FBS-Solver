#pragma once

#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector

namespace FBS {
// a place to aggregate helper functions which do not fit anywhere else
namespace utilities {

	void fillValuesPowerLaw(const double minValue, const double maxValue, std::vector<double>& values, const int power);
	void fillValuesLogarithmic(const double minValue, const double maxValue, std::vector<double>& values);
}
}
