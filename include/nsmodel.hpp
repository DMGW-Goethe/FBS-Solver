#ifndef NSMODEL_HPP
#define NSMODEL_HPP

#include <utility>  // for std::swap

#include "vector.hpp"
#include "eos.hpp"
#include "integrator.hpp"
#include "plotting.hpp"

#define R_INIT 1e-10    // the initial integration radius
#define R_MAX 500.      // the maximum integration radius
#define P_ns_min 1e-15  // the minimum pressure for the "boundary" of the NS

/* this is an abstract class that is supposed to be the backbone for
 * a physical model of a neutron star
 * */
class NSmodel {
public:
    std::shared_ptr<EquationOfState> EOS;   // shared pointer so that multiple NS models can use the same table (saves memory)

    NSmodel(std::shared_ptr<EquationOfState> EOS) : EOS(EOS) {}
    virtual vector dy_dt(const double r, const vector& vars)const  =0;
    static vector dy_dt_static(const double r, const vector& y, const void* params); // this is a static function that receive the class pointer in params and calls dy_dt

    friend std::ostream& operator<<(std::ostream&, const NSmodel&);
    static std::vector<std::string> labels();
};


// this is a class modeling a fermion boson star, akin to
// https://arxiv.org/pdf/2110.11997.pdf
// constructor: EOS (ptr), mu (double), lambda (double), omega (double)
class FermionBosonStar : public NSmodel {
protected:
    vector initial_conditions;
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);

public:
    double mu, lambda, omega;   // holds the defining values of the bosonic scalar field. paricle mass mu, self-interaction parameter lambda, frequency omega
    double rho_0, phi_0;
    // total mass M_T; number of bosons N_B; number of fermions N_F; bosonic radius R_B; fermionic radius R_F; fermionic radius where pressure is zero R_F_0
    double M_T, N_B, N_F, R_B, R_F, R_F_0;

    FermionBosonStar(std::shared_ptr<EquationOfState> EOS, double mu, double lambda=0., double omega=0.)
            : NSmodel(EOS), mu(mu),lambda(lambda), omega(omega), rho_0(0.), phi_0(0.), M_T(0.), N_B(0.), N_F(0.), R_B(0.), R_F(0.), R_F_0(0.) {}

    vector dy_dt(const double r, const vector& vars) const;  // holds the system of ODEs for the Fermion Boson Star
    void set_initial_conditions(const double rho_0, const double phi_0); // holds the FBS init conditions
    virtual int integrate(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), double r_init=R_INIT, double r_end=R_MAX) const ;
    int bisection(double omega_0, double omega_1, int n_mode=0, int max_step=500, double delta_omega=1e-15);
    void shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, double omega_0, double omega_1, int n_mode=0, int max_step=500, double delta_omega=1e-15);
    void evaluate_model(std::vector<integrator::step>& results, integrator::IntegrationOptions intOpts= integrator::IntegrationOptions(), std::string filename="");
    void evaluate_model();

    friend std::ostream& operator<<(std::ostream&, const FermionBosonStar&);
    static std::vector<std::string> labels();

    static const integrator::Event M_converged, Psi_diverging, phi_negative, phi_positive;

};

class FermionBosonStarTLN : public FermionBosonStar {
protected:
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);
    void TLNeq(const double r, const vector& vars, vector& dvars_dr) const;

    FermionBosonStarTLN(const FermionBosonStar& fbs, int H_index, int dH_index, int phi_1_index, int dphi_1_index) :
        FermionBosonStar(fbs), H_index(H_index), dH_index(dH_index), phi_1_index(phi_1_index), dphi_1_index(dphi_1_index), H_0(0.), phi_1_0(0.), lambda_tidal(0.), k2(0.), y_max(0.)
        {
            this->set_initial_conditions();
        }

    const int H_index, dH_index, phi_1_index, dphi_1_index;

public:

    double H_0, phi_1_0;
    double lambda_tidal, k2, y_max;

    FermionBosonStarTLN(const FermionBosonStar& fbs) : FermionBosonStarTLN(fbs, 5, 6, 7, 8) {}

    vector dy_dt(const double r, const vector& vars) const;  // holds the system of ODEs for the Fermion Boson Star + TLN
    virtual void set_initial_conditions(const double phi_1_0=1., const double H_0=1., const double r_init=R_INIT);
    using FermionBosonStar::evaluate_model;
    void evaluate_model(std::vector<integrator::step>& results, integrator::IntegrationOptions intOpts= integrator::IntegrationOptions(), std::string filename="");
    void evaluate_model();
    int bisection_phi_1(double phi_1_0, double phi_1_1, int n_mode=0, int max_step=200, double delta_phi_1=1e-12);

    friend std::ostream& operator<<(std::ostream&, const FermionBosonStarTLN&);
    static std::vector<std::string> labels();

};



class FermionBosonStarTLNInterp : public FermionBosonStarTLN {
protected:
    //void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);
    std::vector<integrator::step> fbs_integration;
    void interpolate(const double r, vector& y_and_dy) const;

public:

    FermionBosonStarTLNInterp(const FermionBosonStar& fbs) : FermionBosonStarTLN(fbs, 0, 1, 2, 3)
    {    integrator::IntegrationOptions intOpts; intOpts.save_intermediate=true; std::vector<integrator::Event> events;
            fbs.integrate(fbs_integration, events, intOpts);   }

    vector dy_dt(const double r, const vector& vars) const;
    void evaluate_model(std::vector<integrator::step>& results, integrator::IntegrationOptions intOpts= integrator::IntegrationOptions(), std::string filename="");
    void set_initial_conditions(const double phi_1_0=1., const double H_0=1., const double r_init=R_INIT);
    int integrate(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), double r_init=R_INIT, double r_end=R_MAX) const ;

};

#endif
