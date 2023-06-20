#ifndef FBS_HPP
#define FBS_HPP


#include "vector.hpp"
#include "eos.hpp"
#include "plotting.hpp"
#include "nsmodel.hpp"
#include "globalvars.hpp"


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

    int find_bosonic_convergence(std::vector<integrator::step>& results, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOps, NUMERIC& R_B_0, bool force=false, NUMERIC r_init=-1., NUMERIC r_end=-1.) const;
    /* Integrates the star until the bosonic field is sufficiently converged phi/phi_0 < PHI_converged, pauses the integration, sets phi=0, and continues the integration
     * to avoid the divergence that is otherwise present due to numerical properties of the system. Only call when omega is found after the bisection! */
    int integrate_and_avoid_phi_divergence(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), bool force = false, std::vector<int> additional_zero_indices={}, NUMERIC r_init=-1., NUMERIC r_end=-1.);

    int bisection_converge_through_infty_behavior(NUMERIC omega_0, NUMERIC omega_1, int n_mode, int max_steps, NUMERIC delta_omega, int verbose);
    int bisection_find_mode(NUMERIC& omega_0, NUMERIC& omega_1, int n_mode, int max_steps, int verbose);
    int bisection_expand_range(NUMERIC& omega_0, NUMERIC& omega_1, int n_mode, int& n_roots_0, int& n_roots_1, int verbose);

    public:
    NUMERIC mu, lambda, omega;
    NUMERIC rho_0, phi_0;
    NUMERIC M_T, N_B, N_F, R_B, R_B_0, R_F, R_F_0, R_G;

    /* Constructor for the FBS class, just sets the relevant values of the class */
    FermionBosonStar(std::shared_ptr<EquationOfState> EOS, NUMERIC mu, NUMERIC lambda=0., NUMERIC omega=0., NUMERIC rho_0=0., NUMERIC phi_0=0.)
            : NSmodel(EOS), mu(mu), lambda(lambda), omega(omega), rho_0(rho_0), phi_0(phi_0), M_T(0.), N_B(0.), N_F(0.), R_B(0.), R_B_0(0.), R_F(0.), R_F_0(0.), R_G(0.) {}

    /* The differential equations describing the FBS. The quantities are a, alpha, P, phi, and Psi, as described in https://arxiv.org/pdf/2110.11997.pdf */
    vector dy_dr(const NUMERIC r, const vector& vars) const;

    /* This function requires mu, lambda, rho_0, phi_0 to be set. It finds the corresponding eigenfrequency omega for the nth mode.
     * omega_0, and omega_1 describe a range in which omega is expected, but the function can extend that range if found to be insufficient*/
    int bisection(NUMERIC omega_0, NUMERIC omega_1, int n_mode=0, int max_step=200, NUMERIC delta_omega=1e-16, int verbose=0);

    /* The initial conditions for a, alpha, phi, Psi, and P */
    virtual vector get_initial_conditions(NUMERIC r_init=-1.) const;

    /* This requires mu, lambda, and rho_0 to be set. It finds phi_0, omega, such that N_B/(N_F+N_B) = NbNf_ratio */
    void shooting_NbNf_ratio(NUMERIC NbNf_ratio, NUMERIC NbNf_accuracy, NUMERIC omega_0, NUMERIC omega_1, int n_mode=0, int max_step=200, NUMERIC delta_omega=1e-16_num, int verbose = 0);

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


#endif
