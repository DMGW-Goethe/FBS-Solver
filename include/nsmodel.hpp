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
    virtual vector dy_dt(const double r, const vector& vars) =0;
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

    vector dy_dt(const double r, const vector& vars);  // holds the system of ODEs for the Fermion Boson Star
    void set_initial_conditions(const double rho_0, const double phi_0); // holds the FBS init conditions
    int integrate(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), double r_init=R_INIT, double r_end=R_MAX) const ;
    int bisection(double omega_0, double omega_1, int n_mode=0, int max_step=500, double delta_omega=1e-15);
    void shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, double omega_0, double omega_1, int n_mode=0, int max_step=500, double delta_omega=1e-15);
    void evaluate_model(std::vector<integrator::step>& results, std::string filename="");
    void evaluate_model();

    friend std::ostream& operator<<(std::ostream&, const FermionBosonStar&);
    static std::vector<std::string> labels();

    static const integrator::Event M_converged, Psi_diverging, phi_negative, phi_positive;

};

class FermionBosonStarTLN : public FermionBosonStar {
protected:
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);

public:

    double H_0, phi_1_0;
    double lambda_tidal, k2, y_max;

    FermionBosonStarTLN(std::shared_ptr<EquationOfState> EOS, double mu, double lambda, double omega)
        : FermionBosonStar(EOS, mu, lambda, omega), H_0(0.), phi_1_0(0.), lambda_tidal(0.), k2(0.), y_max(0.) {}
    FermionBosonStarTLN(const FermionBosonStar& fbs) : FermionBosonStar(fbs) { this->set_initial_conditions(); }

    vector dy_dt(const double r, const vector& vars);  // holds the system of ODEs for the Fermion Boson Star + TLN
    void set_initial_conditions(const double phi_1_0=1., const double H_0=1., const double r_init=R_INIT);
    using FermionBosonStar::evaluate_model;
    void evaluate_model(std::vector<integrator::step>& results, std::string filename="");
    void evaluate_model();
    int bisection_phi_1(double phi_1_0, double phi_1_1, int n_mode=0, int max_step=200, double delta_phi_1=1e-12);

    friend std::ostream& operator<<(std::ostream&, const FermionBosonStarTLN&);
    static std::vector<std::string> labels();

    static const integrator::Event dphi_1_diverging, phi_1_negative, phi_1_positive;
};


// this is a class modeling a fermion boson star in the two-fluid ansatz, akin to
// PHYSICAL REVIEW D 105, 123010 (2022)
// generally, this "two fluid FBS" can describe any neutron star made out of two non-coupling fluids (only interact gravitationally) with each their EoS.
// it must not necessarily be neuron-matter + dark matter, but can be a combination of two arbitrary fluids with EoS.
// constructor: EOS (ptr), EOS2 (ptr), mu (double), lambda (double)
class TwoFluidFBS : public NSmodel {
protected:
    vector initial_conditions;
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);

public:
    double mu, lambda;   // holds the defining values of the bosonic scalar field. paricle mass mu, self-interaction parameter lambda
    double rho1_0, rho2_0;	// initial conditions, central density of fluid 1 and 2 respectively
    // total mass M_T; total mass of fluid 1 (2): M_1  (M_2); radius R_1 , R_2 (99% of matter included); radius where pressure is zero R_1_0 , R_2_0;
    double M_T, M_1, M_2, R_1, R_1_0, R_2, R_2_0, k2, lambda_tidal;

	std::shared_ptr<EquationOfState> EOS_fluid2;	// EOS of the second fluid

    TwoFluidFBS(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2)
            : NSmodel(EOS1), EOS_fluid2(EOS2), mu(mu),lambda(lambda), rho1_0(0.), rho2_0(0.), M_T(0.), M_1(0.), M_2(0.), R_1(0.), R_1_0(0.), R_2(0.), R_2_0(0.), k2(0.), lambda_tidal(0.) {}

    vector dy_dt(const double r, const vector& vars);  // holds the system of ODEs for the Fermion Boson Star
    void set_initial_conditions(const double rho1_0, const double rho2_0); // holds the FBS init conditions
    int integrate(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), double r_init=R_INIT, double r_end=R_MAX) const ;
    // void shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, int max_step=500);
    void evaluate_model(std::vector<integrator::step>& results, std::string filename="");
    void evaluate_model();

    friend std::ostream& operator<<(std::ostream&, const TwoFluidFBS&);
    static std::vector<std::string> labels();

    static const integrator::Event all_Pressure_zero;
};

#endif
