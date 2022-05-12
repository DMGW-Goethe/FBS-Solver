#ifndef NSMODEL_HPP
#define NSMODEL_HPP

#include "vector.hpp"
#include "eos.hpp"


/* this is an abstract class that is supposed to be the backbone for
 * a physical model of a neutron star
 * */
class NSmodel {
public:
    std::shared_ptr<EquationOfState> EOS;

    NSmodel(std::shared_ptr<EquationOfState> EOS) : EOS(EOS) {}
    virtual vector dy_dt(const double r, const vector& vars) =0;
    static vector dy_dt_static(const double r, const vector& y, const void* params); // this is a static function that receive the class pointer in params and calls dy_dt
};


/* this is a class modeling a fermion boson star, akin to
 * https://arxiv.org/pdf/2110.11997.pdf
 * */
class FermionBosonStar : public NSmodel {
public:
    double mu, lambda, omega;

    FermionBosonStar(std::shared_ptr<EquationOfState> EOS, double mu, double lambda, double omega) : NSmodel(EOS), mu(mu),lambda(lambda), omega(omega) {}

    vector dy_dt(const double r, const vector& vars) ;
    vector initial_conditions(const double r0, const double rho_0, const double phi_0) ;
};




#endif
