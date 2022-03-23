#ifndef NSMODEL_HPP
#define NSMODEL_HPP

#include "vector.hpp"
#include "eos.hpp"



class NSmodel {

public:
    double mu, lambda, omega;
    std::shared_ptr<EquationOfState> EOS;

    NSmodel(std::shared_ptr<EquationOfState> EOS, double mu, double lambda, double omega) : mu(mu),lambda(lambda), omega(omega), EOS(EOS) {}

    vector dy_dt(const double r, const vector& vars) ;
    static vector dy_dt_static(const double r, const vector& y, const void* params); // this is a static function that receive the class pointer in params and calls dy_dt

};



#endif
