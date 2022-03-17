#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;


class vector : public ublas::vector<double> {
public:
    using ublas::vector<double>::vector;

    vector(std::initializer_list<double> list);
};
#endif
