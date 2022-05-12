#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/io.hpp"

namespace ublas = boost::numeric::ublas;

/* this class is simply a wrapper for the ublas::vector class
 * that includes a convenient constructor and some nan checks */
class vector : public ublas::vector<double> {
public:
    using ublas::vector<double>::vector;

    vector(std::initializer_list<double> list);

    static bool is_nan(const vector& v);
};

#endif
