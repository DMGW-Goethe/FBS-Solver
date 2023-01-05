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
#define PHI_converged 1e-4
#define INT_converged 1e-7
#define M_T_converged 1e-15

/*  NSmodel
 * this is an abstract class that is supposed to be the backbone for
 *  a physical model of a neutron star
 * This implementation allows the integration of the differential equations
 *  describing the neutron star
 * */
class NSmodel {
protected:
    double r_init, r_end;
public:
    /* A pointer to the Equation of State describing the neutron star matter */
    std::shared_ptr<EquationOfState> EOS;

    /* Constructor */
    NSmodel(std::shared_ptr<EquationOfState> EOS) : r_init(R_INIT), r_end(R_MAX), EOS(EOS) {}

    /* This function gives the derivatives for the differential equations */
    virtual vector dy_dr(const double r, const vector& vars) const = 0;

    /* This is a static wrapper function, that calls the dy_dr function */
    static vector dy_dr_static(const double r, const vector& y, const void* params);
    /*static vector dy_dr_static_scaled(const double r, const vector& y, const void* params);*/

    /* The initial conditions for a, alpha, phi, Psi, and P */
    virtual vector get_initial_conditions(double r_init=R_INIT) const = 0;

    /* This function calls the integrator and returns the results of the integration */
    int integrate(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, const vector initial_conditions, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), double r_init=-1., double r_end=-1.) const;
    /* int scaled_integration(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, const vector initial_conditions, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), double r_init=R_INIT, double r_end=R_MAX) const;*/

    /* For easy output the class should define an << operator */
    friend std::ostream& operator<<(std::ostream&, const NSmodel&);
    /* The labels for the different parameters that are output by the << operator */
    static std::vector<std::string> labels();
};


/* FermionBosonStar
 * This class models a fermion boson star (FBS)
 *  akin to https://arxiv.org/pdf/2110.11997.pdf
 * The star consists of neutron matter and a bosonic field phi
 *  and is characterized by the parameters
 *      mu      : mass of the boson field
 *      lambda  : strength of the quartic field self-interaction
 *      omega   : the eigenfrequency of the field
 *      rho_0   : The central density of the neutron matter
 *      phi_0   : The central field value of the phi field
 *
 * The values mu, lambda, omega, phi_0 need to be related such that phi->0 at infty
 * To obtain omega for given mu, lambda, phi_0, the function bisection finds omega to an approximate degree
 * Alternatively NbNf_ratio can be used to obtain the pair of phi_0, omega such that N_B/N_F is some specific value
 *
 * Once all parameters have been set, the evaluate_model function calculates the properties of the FBS
 *      M_T     : The ADM mass
 *      N_B     : The particle number of the bosonic field
 *      N_F     : The particle number of the neutron star matter
 *      R_B     : The bosonic radius, at which > 0.99N_B field particles are at
 *      R_F_0   : The fermionic radius, where P < P_ns_min
 */
class FermionBosonStar : public NSmodel {
protected:
    /* Calculates the parameters M_T, N_B, N_F, R_B, R_F_0 for a given integration contained in results,events */
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);

    int find_bosonic_convergence(std::vector<integrator::step>& results, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOps, double& R_B_0, bool force=false, double r_init=-1., double r_end=-1.) const;
    /* Integrates the star until the bosonic field is sufficiently converged phi/phi_0 < PHI_converged, pauses the integration, sets phi=0, and continues the integration
     * to avoid the divergence that is otherwise present due to numerical properties of the system. Only call when omega is found after the bisection! */
    int integrate_and_avoid_phi_divergence(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), bool force = false, std::vector<int> additional_zero_indices={}, double r_init=-1., double r_end=-1.);

    public:
    double mu, lambda, omega;
    double rho_0, phi_0;
    double M_T, N_B, N_F, R_B, R_B_0, R_F, R_F_0;

    /* Constructor for the FBS class, just sets the relevant values of the class */
    FermionBosonStar(std::shared_ptr<EquationOfState> EOS, double mu, double lambda=0., double omega=0., double rho_0=0., double phi_0=0.)
            : NSmodel(EOS), mu(mu), lambda(lambda), omega(omega), rho_0(rho_0), phi_0(phi_0), M_T(0.), N_B(0.), N_F(0.), R_B(0.), R_B_0(0.), R_F(0.), R_F_0(0.) {}

    /* The differential equations describing the FBS. The quantities are a, alpha, P, phi, and Psi, as described in https://arxiv.org/pdf/2110.11997.pdf */
    vector dy_dr(const double r, const vector& vars) const;

    /* This function requires mu, lambda, rho_0, phi_0 to be set. It finds the corresponding eigenfrequency omega for the nth mode.
     * omega_0, and omega_1 describe a range in which omega is expected, but the function can extend that range if found to be insufficient*/
    int bisection(double omega_0, double omega_1, int n_mode=0, int max_step=200, double delta_omega=1e-16, int verbose=0);

    /* The initial conditions for a, alpha, phi, Psi, and P */
    virtual vector get_initial_conditions(double r_init=-1.) const;

    /* This requires mu, lambda, and rho_0 to be set. It finds phi_0, omega, such that N_B/N_F = NbNf_ratio */
    void shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, double omega_0, double omega_1, int n_mode=0, int max_step=200, double delta_omega=1e-16);

    /* Integrates the DE while avoiding the phi divergence and calculates the FBS properties
     * Returns the results and optionally ouputs them into a file*/
    void evaluate_model(std::vector<integrator::step>& results, integrator::IntegrationOptions intOpts= integrator::IntegrationOptions(), std::string filename="");
    /* Wrapper if the results are not wanted */
    void evaluate_model();

    /* This function outputs parameters and properties, in the order given by the labels function */
    friend std::ostream& operator<<(std::ostream&, const FermionBosonStar&);
    static std::vector<std::string> labels();

    /* These events are used for different integration purposes */
    static const integrator::Event M_converged, Psi_diverging, phi_negative, phi_positive, phi_converged, integration_converged, P_min_reached, Psi_positive;

};


/* FermionBosonStarTLN
 * This class models a perturbed fermion boson star (FBS). Here, the perturbations
 *  are described by the metric perturbation H(r) and bosonic field perturbation phi_1(r)
 *
 * The class inherits all parameters mu, lambda, rho_0, phi_0, omega
 *  and properties M_T, N_B, N_F, R_B, R_F_0 from the parent class
 *
 *  The additional parameters are
 *      H_0     : The central value of the metric perturbations
 *      phi_1_0 : The central value of the field perturbation
 *
 *  The additional properties are
 *      lambda_tidal  : The tidal deformability as given by eq (26) of https://arxiv.org/pdf/1704.08651.pdf
 *      k2            : The tidal love number as given by eq (23) of https://arxiv.org/pdf/0711.2420.pdf
 *      y_max         : The value of y = r H'/H at the maximum that goes into the above equations
 *      R_ext         : The point at which this maximum was found and the values computed
 *
 *  The implementation assumes that the FermionBosonStar that is used to construct this class already has
 *   all parameters consistent and properties calculated!
 *
 *  The normalization of H, phi_1 is arbitrary and only affects each other. Therefore, H_0 = 1 usually
 *  The corresponding phi_1_0 can be found with another bisection algorithm
 */
class FermionBosonStarTLN : public FermionBosonStar {
protected:
    /* Calculates the parameters lambda_tidal, k2, y_max, R_ext for a given integration contained in results,events */
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);

    /* Integrates the star until the bosonic field is sufficiently converged phi/phi_0 < PHI_converged, pauses the integration, sets phi=0, phi_1 = 0, and continues the integration
     * to avoid the divergence that is otherwise present due to numerical properties of the system. For phi=0 the DE for H,phi_1 decouple and we are only interested in H anyways. */
    /*int integrate_and_avoid_phi_divergence(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), double r_init=R_INIT, double r_end=R_MAX) const;*/

public:
    double H_0, phi_1_0;
    double lambda_tidal, k2, y_max, R_ext;

    /*FermionBosonStarTLN(std::shared_ptr<EquationOfState> EOS, double mu, double lambda, double omega)
        : FermionBosonStar(EOS, mu, lambda, omega), H_0(1.), phi_1_0(0.), lambda_tidal(0.), k2(0.), y_max(0.), R_ext(0) {}*/ // TODO: Check, this shouldn't work without call to evaluate_model of the parent class
    /* The constructor for the class. Assumes that the FermionBosonStar has its properties calculated! */
    FermionBosonStarTLN(const FermionBosonStar& fbs) : FermionBosonStar(fbs), H_0(1.), phi_1_0(0.), lambda_tidal(0.), k2(0.), y_max(0.), R_ext(0.) {  }

    /* The differential equations describing the FBS + TLN. The quantities are a, alpha, phi, Psi, P, H, dH, phi_1, dphi_1 */
    vector dy_dr(const double r, const vector& vars) const;

    /* The initial conditions for a, alpha, phi, Psi, P, H, dH, phi_1, dphi_1*/
    vector get_initial_conditions(double r_init=R_INIT) const;

    /* This function requires the FBS parameters and H_0 to be set. It finds the corresponding phi_1_0 */
    int bisection_phi_1(double phi_1_0, double phi_1_1, int n_mode=0, int max_step=200, double delta_phi_1=1e-12, int verbose = 0);

    using FermionBosonStar::evaluate_model; // TODO: Check Necessity

    /* Integrates the DE while avoiding the phi divergence and calculates the FBS properties
     * Returns the results and optionally ouputs them into a file*/
    void evaluate_model(std::vector<integrator::step>& results, std::string filename="");
    void evaluate_model();

    /* This function outputs parameters and properties, in the order given by the labels function */
    friend std::ostream& operator<<(std::ostream&, const FermionBosonStarTLN&);
    static std::vector<std::string> labels();

    /* These events are used for different integration purposes */
    static const integrator::Event dphi_1_diverging, phi_1_negative, phi_1_positive;
};


// this is a class modeling a fermion boson star in the two-fluid ansatz, akin to
// PHYSICAL REVIEW D 105, 123010 (2022)
// generally, this "two fluid FBS" can describe any neutron star made out of two non-coupling fluids (only interact gravitationally) with each their EoS.
// it must not necessarily be neuron-matter + dark matter, but can be a combination of two arbitrary fluids with EoS.
// constructor: EOS (ptr), EOS2 (ptr), mu (double), lambda (double)
class TwoFluidFBS : public NSmodel {
protected:
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);

public:
    //double mu, lambda;   // holds the defining values of the bosonic scalar field. paricle mass mu, self-interaction parameter lambda
    double rho1_0, rho2_0;	// initial conditions, central density of fluid 1 and 2 respectively
    // total mass M_T; total mass of fluid 1 (2): M_1  (M_2); radius R_1 , R_2 (99% of matter included); radius where pressure is zero R_1_0 , R_2_0;
    double M_T, M_1, M_2, R_1, R_1_0, R_2, R_2_0, C, k2, lambda_tidal;
	double N_1, N_2; // fermion/boson numbers computed using the noether current

	// define variables in case we use the effective bosonic EoS:
	bool use_effective_bosonicEoS = false;
	double mu = 1., lambda = 1.;

	std::shared_ptr<EquationOfState> EOS_fluid2;	// EOS of the second fluid

    TwoFluidFBS()
            : NSmodel(nullptr), EOS_fluid2(nullptr), rho1_0(0.), rho2_0(0.), M_T(0.), M_1(0.), M_2(0.), R_1(0.), R_1_0(0.), R_2(0.), R_2_0(0.), k2(0.), lambda_tidal(0.), N_1(0.), N_2(0.), mu(0.), lambda(0.), use_effective_bosonicEoS(false) {}

	TwoFluidFBS(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2)
        : NSmodel(EOS1), EOS_fluid2(EOS2), rho1_0(0.), rho2_0(0.), M_T(0.), M_1(0.), M_2(0.), R_1(0.), R_1_0(0.), R_2(0.), R_2_0(0.), k2(0.), lambda_tidal(0.), N_1(0.), N_2(0.) {}

    TwoFluidFBS(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2, double mu, double lambda)
            : NSmodel(EOS1), EOS_fluid2(EOS2), rho1_0(0.), rho2_0(0.), M_T(0.), M_1(0.), M_2(0.), R_1(0.), R_1_0(0.), R_2(0.), R_2_0(0.), k2(0.), lambda_tidal(0.), N_1(0.), N_2(0.), mu(mu), lambda(lambda), use_effective_bosonicEoS(true) {}

    //vector dy_dr(const double r, const vector& vars);  // holds the system of ODEs for the Fermion Boson Star
	/* The differential equations describing the two-fluid FBS. The quantities are nu, m1, m2, P1, P2 and y as described in PHYS. REV. D 105, 123010 (2022) */
    vector dy_dr(const double r, const vector& vars) const;

    vector get_initial_conditions(const double r_init=R_INIT) const; // holds the FBS init conditions
   // void shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, int max_step=500);
    void evaluate_model(std::vector<integrator::step>& results, std::string filename="");
    void evaluate_model();

    friend std::ostream& operator<<(std::ostream&, const TwoFluidFBS&);
    static std::vector<std::string> labels();

    static const integrator::Event all_Pressure_zero;
};

#endif
