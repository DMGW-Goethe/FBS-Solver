#include "vector.hpp"



vector::vector(std::initializer_list<double> list) : vector(list.size())  {
        for(unsigned i = 0; i < list.size(); ++ i)
            (*this)[i] = *(list.begin() + i);
    }
