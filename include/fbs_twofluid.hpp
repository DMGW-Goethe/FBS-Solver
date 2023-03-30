#ifndef FBS_TWOFLUID_HPP
#define FBS_TWOFLUID_HPP

#include <utility>  // for std::swap

#include "vector.hpp"
#include "eos.hpp"
#include "integrator.hpp"
#include "ns_twofluid.hpp"
#include "nsmodel.hpp"


// this is a class modeling a fermion boson star in the two-fluid ansatz, akin to
// PHYSICAL REVIEW D 105, 123010 (2022)
// generally, this "two fluid FBS" can describe any neutron star made out of two non-coupling fluids (only interact gravitationally) with each their EoS.
// it must not necessarily be neuron-matter + dark matter, but can be a combination of two arbitrary fluids with EoS.
// constructor: EOS (ptr), EOS2 (ptr), mu (double), lambda (double)
class TwoFluidFBS : public NSTwoFluid {

public:
	// define variables in case we use the effective bosonic EoS:
	double mu = 1., lambda = 1.;

	//std::shared_ptr<EquationOfState> EOS_fluid2;	// EOS of the second fluid

    TwoFluidFBS()
            : NSTwoFluid(nullptr, nullptr), mu(1.), lambda(1.) {}

    TwoFluidFBS(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2, double mu, double lambda)
            : NSTwoFluid(EOS1, EOS2), mu(mu), lambda(lambda) {}

	/* EOM are identical for any two fluid system, no matter what kind of EoS is used, threfore it is not necessary to re-declare the ODEs here.
	   See the parent class for details about the implementation. */

    vector get_initial_conditions(const double r_init=R_INIT) const; // holds the FBS init conditions

    static std::vector<std::string> labels();
};


#endif