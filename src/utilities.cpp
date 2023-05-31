#include "utilities.hpp"


void utilities::fillValuesPowerLaw(const NUMERIC minValue, const NUMERIC maxValue, std::vector<NUMERIC>& values, const int power)
{
    if(values.size() == 1) {
		values[0] = minValue;
		return;
	}
	
	if(power == 1)
    {
        const NUMERIC dValue = (NUMERIC)(maxValue - minValue) / (NUMERIC)(values.size() - 1);
        for(size_t i = 0; i < values.size(); i++)
            values[i] = minValue + dValue * i;

        return;
    }

    fillValuesPowerLaw(0._num, 1._num, values, 1);

    for(size_t i = 0; i < values.size(); i++)
    {
        values[i] = pow(values[i], power);
        values[i] *= maxValue - minValue;
        values[i] += minValue;
    }
}

void utilities::fillValuesLogarithmic(const NUMERIC minValue, const NUMERIC maxValue, std::vector<NUMERIC>& values)
{
    fillValuesPowerLaw(log(minValue), log(maxValue), values, 1);

    for(size_t i = 0; i < values.size(); i++)
        values[i] = exp(values[i]);
}