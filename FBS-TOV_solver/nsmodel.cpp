#include "nsmodel.hpp"

vector NSmodel::dy_dt_static(const double r, const vector& y, const void* params) {
    NSmodel* m = (NSmodel*)params;
    return m->dy_dt(r, y);
}

// define the system of coupled ODEs for a Fermion-Boson Star:
vector FermionBosonStar::dy_dt(const double r, const vector& vars) {

    // rename input & class variables for simpler use:
    const double a = vars[0]; const double alpha = vars[1]; const double Phi = vars[2]; const double Psi = vars[3]; double P = vars[4];
    EquationOfState& myEOS = *(this->EOS);
    const double mu = this->mu; const double lambda = this->lambda; const double omega = this->omega;

    // define hydrodynamic quantities:
    double rho = 1.;      // restmass density, must be set using EOS
    double epsilon = 1.;  // specific energy denstiy, must be set either through EOS or hydrodynamic relations
    // epsilon is related to the total energy density "e" by: e = rho*(1+epsilon)

    if(P < 0.) P = 0.;  // need this to prevent NaN errors...

    // apply the EOS:
    myEOS.callEOS(rho, epsilon, P); // change rho and epsilon by pointer using EOS member function

    if( vector::is_nan(vars)) {
        std::cout << "Nan found: " << vars << std::endl;
        assert(false);}

    // compute the ODEs:
    double da_dr = 0.5* a * ( (1.-a*a) / r + 4.*M_PI*r*( (omega*omega/ alpha/alpha + mu*mu + 0.5*lambda*Phi*Phi )*a*a*Phi*Phi + Psi*Psi + 2.*a*a*rho*(1.+epsilon) ) );
    double dalpha_dr = 0.5* alpha * ( (a*a-1.) / r + 4.*M_PI*r*( (omega*omega/ alpha/alpha - mu*mu - 0.5*lambda*Phi*Phi )*a*a*Phi*Phi + Psi*Psi + 2.*a*a*P ) );
    double dPhi_dr = Psi;
    double dPsi_dr = -( 1. + a*a - 4.*M_PI*r*r*a*a*( mu*mu*Phi*Phi + 0.5*lambda*Phi*Phi*Phi*Phi + rho*(1.+epsilon) - P ))*Psi/r - (omega*omega/ alpha/alpha - mu*mu - lambda*Phi*Phi )*a*a*Phi;
    double dP_dr = -(rho*(1.+epsilon) + P)*dalpha_dr/alpha;

    // write the ODE values into output vector:
    return vector({da_dr, dalpha_dr, dPhi_dr, dPsi_dr, dP_dr});
}


vector FermionBosonStar::initial_conditions(const double r0, const double rho_0, const double phi_0) {
        return vector( {1.0, 1.0, phi_0, 0., this->EOS->get_P_from_rho(rho_0)});
}

