#ifndef FBS_TLN_HPP
#define FBS_TLN_HPP

#include <utility>  // for std::swap

#include "vector.hpp"
#include "eos.hpp"
#include "plotting.hpp"
#include "nsmodel.hpp"
#include "fbs.hpp"

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


    int bisection_phi_1_find_mode(NUMERIC& phi_1_0_l, NUMERIC& phi_1_0_r, int n_mode, int max_steps, int verbose);
    int bisection_phi_1_converge_through_infty_behavior(NUMERIC& phi_1_0_l, NUMERIC& phi_1_0_r, int max_steps, NUMERIC delta_phi_1, int verbose);


public:
    NUMERIC H_0, phi_1_0;
    NUMERIC lambda_tidal, k2, y_max, R_ext;

    /*FermionBosonStarTLN(std::shared_ptr<EquationOfState> EOS, NUMERIC mu, NUMERIC lambda, NUMERIC omega)
        : FermionBosonStar(EOS, mu, lambda, omega), H_0(1.), phi_1_0(0.), lambda_tidal(0.), k2(0.), y_max(0.), R_ext(0) {}*/ // TODO: Check, this shouldn't work without call to evaluate_model of the parent class
    /* The constructor for the class. Assumes that the FermionBosonStar has its properties calculated! */
    FermionBosonStarTLN(const FermionBosonStar& fbs) : FermionBosonStar(fbs), H_0(1.), phi_1_0(0.), lambda_tidal(0.), k2(0.), y_max(0.), R_ext(0.) {  }

    /* The differential equations describing the FBS + TLN. The quantities are a, alpha, phi, Psi, P, H, dH, phi_1, dphi_1 */
    vector dy_dr(const NUMERIC r, const vector& vars) const;

    /* The initial conditions for a, alpha, phi, Psi, P, H, dH, phi_1, dphi_1*/
    vector get_initial_conditions(NUMERIC r_init=R_INIT) const;

    /* This function requires the FBS parameters and H_0 to be set. It finds the corresponding phi_1_0 */
    int bisection_phi_1(NUMERIC phi_1_0, NUMERIC phi_1_1, int n_mode=0, int max_step=200, NUMERIC delta_phi_1=1e-12_num, int verbose = 0);

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


#endif
