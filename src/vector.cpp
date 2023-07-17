#include "vector.hpp"

using namespace FBS;

// constructor
vector::vector(std::initializer_list<double> list) : vector(list.size())  {
        for(unsigned i = 0; i < list.size(); ++ i)
            (*this)[i] = *(list.begin() + i);
}

vector vector::sub_range(int begin, int end) const  {
    //vector s(this->begin()+begin, this->begin() + begin + end);
    //ublas::vector_range<vector> vr((const vector&)*this, ublas::range(begin, end));
    //return vector(vr);
    vector vr(end-begin);
    for(int i = begin; i < end; i++)
        vr[i-begin] = (*this)[i];
    return vr;
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
