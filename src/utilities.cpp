#include "utilities.hpp"

using namespace FBS;

void utilities::fillValuesPowerLaw(const double minValue, const double maxValue, std::vector<double>& values, const int power)
{
    if(values.size() == 1) {
		values[0] = minValue;
		return;
	}

	if(power == 1)
    {
        const double dValue = double(maxValue - minValue) / double(values.size() - 1);
        for(size_t i = 0; i < values.size(); i++)
            values[i] = minValue + dValue * i;

        return;
    }

    fillValuesPowerLaw(0.0, 1.0, values, 1);

    for(size_t i = 0; i < values.size(); i++)
    {
        values[i] = pow(values[i], power);
        values[i] *= maxValue - minValue;
        values[i] += minValue;
    }
}

void utilities::fillValuesLogarithmic(const double minValue, const double maxValue, std::vector<double>& values)
{
    fillValuesPowerLaw(log(minValue), log(maxValue), values, 1);

    for(size_t i = 0; i < values.size(); i++)
        values[i] = exp(values[i]);
}
