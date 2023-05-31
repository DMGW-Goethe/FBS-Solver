#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
#include "vector.hpp"

// a place to aggregate helper functions which do not fit anywhere else
namespace utilities {

	void fillValuesPowerLaw(const NUMERIC minValue, const NUMERIC maxValue, std::vector<NUMERIC>& values, const int power);
	void fillValuesLogarithmic(const NUMERIC minValue, const NUMERIC maxValue, std::vector<NUMERIC>& values);
}

#endif