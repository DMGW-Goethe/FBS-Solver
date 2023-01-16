#ifndef VECTOR_HPP
#define VECTOR_HPP

// include the external boost/ublas library and use the n-dimensional vector class
// see: https://www.boost.org/
#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
//#include <boost/range.hpp>

namespace ublas = boost::numeric::ublas;

/* this class is simply a wrapper for the ublas::vector class
 * that includes a convenient constructor and some nan checks */
class vector : public ublas::vector<double> {
public:
    using ublas::vector<double>::vector;

    vector(std::initializer_list<double> list); // created a vector with size len(list)

    vector sub_range(int begin, int end) const ;

    static bool is_nan(const vector& v);    // functions to check for NaNs
};

#endif
