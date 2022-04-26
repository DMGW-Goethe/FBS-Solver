#include "vector.hpp"



vector::vector(std::initializer_list<double> list) : vector(list.size())  {
        for(unsigned i = 0; i < list.size(); ++ i)
            (*this)[i] = *(list.begin() + i);
    }

bool vector::is_nan(const vector& v) {
    for(unsigned i =0; i < v.size(); i++) {
        if(std::isnan(v[i]))
            return true;
    }
    return false;
}

