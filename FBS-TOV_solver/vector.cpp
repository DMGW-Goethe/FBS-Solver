#include "vector.hpp"

// constructor
vector::vector(std::initializer_list<double> list) : vector(list.size())  {
        for(unsigned i = 0; i < list.size(); ++ i)
            (*this)[i] = *(list.begin() + i);
}

// NaN-checker function
// returns NaN if any entry is NaN
bool vector::is_nan(const vector& v) {
    for(unsigned i =0; i < v.size(); i++) {
        if(std::isnan(v[i]))
            return true;
    }
    return false;
}