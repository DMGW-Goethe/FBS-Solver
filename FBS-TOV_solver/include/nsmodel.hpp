#ifndef NSMODEL_HPP
#define NSMODEL_HPP

#include "vector.hpp"
#include "eos.hpp"
#include "integrator.hpp"
#include "plotting.hpp"


/* this is an abstract class that is supposed to be the backbone for
 * a physical model of a neutron star
 * */
class NSmodel {
public:
    std::shared_ptr<EquationOfState> EOS;   // shared pointer so that multiple NS models can use the same table (saves memory)

    NSmodel(std::shared_ptr<EquationOfState> EOS) : EOS(EOS) {}
    virtual vector dy_dt(const double r, const vector& vars) =0;
    static vector dy_dt_static(const double r, const vector& y, const void* params); // this is a static function that receive the class pointer in params and calls dy_dt
};


// this is a class modeling a fermion boson star, akin to
// https://arxiv.org/pdf/2110.11997.pdf
// constructor: EOS (ptr), mu (double), lambda (double), omega (double)
class FermionBosonStar : public NSmodel {
public:
    double mu, lambda, omega;   // holds the defining values of the bosonic scalar field. paricle mass mu, self-interaction parameter lambda, frequency omega
    vector initial_conditions;
    // total mass M_T; number of bosons N_B; number of fermions N_F; bosonic radius R_B; fermionic radius R_F; fermionic radius where pressure is zero R_F_0
    double M_T, N_B, N_F, R_B, R_F, R_F_0;

    FermionBosonStar(std::shared_ptr<EquationOfState> EOS, double mu, double lambda, double omega) : NSmodel(EOS), mu(mu),lambda(lambda), omega(omega) {}

    vector dy_dt(const double r, const vector& vars);  // holds the system of ODEs for the Fermion Boson Star
    void set_initial_conditions(const double r0, const double rho_0, const double phi_0); // holds the FBS init conditions
    void bisection(double omega_0, double omega_1, int n_mode=0, int max_step=1000, double delta_omega=1e-15);
    void evaluate_model(std::string filename="");
};

#endif
