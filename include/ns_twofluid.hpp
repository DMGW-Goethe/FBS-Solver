#ifndef NS_TWOFLUID_HPP
#define NS_TWOFLUID_HPP

#include <utility>  // for std::swap

#include "vector.hpp"
#include "eos.hpp"
#include "integrator.hpp"
#include "nsmodel.hpp"


// this is a class modeling a fermion boson star in the two-fluid ansatz, akin to
// PHYSICAL REVIEW D 105, 123010 (2022)
// generally, this "two fluid FBS" can describe any neutron star made out of two non-coupling fluids (only interact gravitationally) with each their EoS.
// it must not necessarily be neuron-matter + dark matter, but can be a combination of two arbitrary fluids with EoS.
// constructor: EOS (ptr), EOS2 (ptr)
class NSTwoFluid : public NSmodel {
protected:
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);

public:
	std::shared_ptr<EquationOfState> EOS_fluid2;	// EOS of the second fluid

    //NUMERIC mu, lambda;   // holds the defining values of the bosonic scalar field. paricle mass mu, self-interaction parameter lambda
    NUMERIC rho1_0, rho2_0;	// initial conditions, central density of fluid 1 and 2 respectively
    // total mass M_T; total mass of fluid 1 (2): M_1  (M_2); radius R_1 , R_2 (99% of matter included); radius where pressure is zero R_1_0 , R_2_0;
    NUMERIC M_T, M_1, M_2, R_1, R_1_0, R_2, R_2_0, C, k2, lambda_tidal;
	NUMERIC N_1, N_2; // fermion/boson numbers computed using the noether current



	NSTwoFluid(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2)
        : NSmodel(EOS1), EOS_fluid2(EOS2), rho1_0(0.), rho2_0(0.), M_T(0.), M_1(0.), M_2(0.), R_1(0.), R_1_0(0.), R_2(0.), R_2_0(0.), k2(0.), lambda_tidal(0.), N_1(0.), N_2(0.) {}

    //vector dy_dr(const NUMERIC r, const vector& vars);  // holds the system of ODEs for the Fermion Boson Star
	/* The differential equations describing the two-fluid FBS. The quantities are nu, m1, m2, P1, P2 and y as described in PHYS. REV. D 105, 123010 (2022) */
    vector dy_dr(const NUMERIC r, const vector& vars) const;

    vector get_initial_conditions(const NUMERIC r_init=R_INIT) const; // holds the FBS init conditions
    void evaluate_model(std::vector<integrator::step>& results, std::string filename="");
    void evaluate_model();

    friend std::ostream& operator<<(std::ostream&, const NSTwoFluid&);
    static std::vector<std::string> labels();

    static const integrator::Event all_Pressure_zero;
};


#endif
