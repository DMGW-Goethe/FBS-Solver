#include "fbs_twofluid.hpp"

using namespace FBS;

/***********************
 * TwoFluidFBS *
 ***********************/

// initial conditions for one NS-matter fluid and one effective bosonic fluid
vector TwoFluidFBS::get_initial_conditions(const double r_init) const {
	// effective central energy density of 2nd fluid (uses the effective bosonic EOS).
	// rho2_0 corresponds to the central value of the scalar field phi_0
    double e_phi_c = 2.* std::pow(mu ,2)* std::pow(rho2_0,2) + 1.5*lambda*std::pow(rho2_0,4);
    // d/dr (nu, m1, m2, p1, p2, y(r))
    return vector({0.0, 0.0, 0.0, rho1_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(rho1_0, 0.) : 0.,
                                       e_phi_c > this->EOS_fluid2->min_e() ? this->EOS_fluid2->get_P_from_e(e_phi_c) : 0., 2.0});
}

std::vector<std::string> TwoFluidFBS::labels() {
    // labels for the effective bosonic EOS case:
    return std::vector<std::string>({"M_T", "rho_0", "phi_0", "R_F", "R_F_0", "M_F", "N_F", "R_B", "R_B_0", "M_B", "N_B", "M_B/M_F", "N_B/N_F", "C", "k2", "lambda_tidal"});
}
