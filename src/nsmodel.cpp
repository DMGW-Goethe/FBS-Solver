#include "nsmodel.hpp"

using std::exp;
using std::abs;

vector NSmodel::dy_dt_static(const double r, const vector& y, const void* params) {
    NSmodel* m = (NSmodel*)params;
    return m->dy_dt(r, y);
}

// define the system of coupled ODEs for a Fermion-Boson Star:
vector FermionBosonStar::dy_dt(const double r, const vector& vars) {

    // rename input & class variables for simpler use:
    const double v = vars[0]; const double u = vars[1]; const double phi = vars[2]; const double Psi = vars[3]; double P = vars[4];
    EquationOfState& myEOS = *(this->EOS);
    const double mu = this->mu; const double lambda = this->lambda; const double omega = this->omega;

    // define hydrodynamic quantities:
    double rho = 1.;      // restmass density, must be set using EOS
    double epsilon = 1.;  // specific energy denstiy, must be set either through EOS or hydrodynamic relations
    // epsilon is related to the total energy density "e" by: e = rho*(1+epsilon)
    const double V = mu*mu*phi*phi + lambda/2.*pow(phi, 4);
    const double dV_deps = mu*mu + lambda*phi*phi;
    //const double ddV_deps2 = lambda;

    // apply the EOS:
    if(P <= 0. || P < myEOS.min_P() || P < P_ns_min)  {
        P = 0.; rho = 0.; epsilon = 0.;
    } else {
        myEOS.callEOS(rho, epsilon, P); // change rho and epsilon by reference using EOS member function
    }

    // compute the ODEs:
    const double dv_dr = (exp(u) - 1.)/r + 8.*M_PI*r*( omega*omega*exp(u-v)*phi*phi - V*exp(u) + Psi*Psi + P*exp(u));
    const double du_dr = -(exp(u) - 1.)/r + 8.*M_PI*r*( omega*omega*exp(u-v)*phi*phi + V*exp(u) + Psi*Psi + rho*(1.+epsilon)*exp(u));
    const double dphi_dr = Psi;
    const double dPsi_dr = (dV_deps*exp(u) - omega*omega*exp(u-v))*phi + (du_dr/2. - dv_dr/2. - 2./r)*Psi;
    const double dP_dr = - (P + rho*(1.+epsilon)) * dv_dr/2.;

    // write the ODE values into output vector:
    return vector({dv_dr, du_dr, dphi_dr, dPsi_dr, dP_dr});
}


const integrator::Event FermionBosonStar::M_converged = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void *params) {
                                                                                                        const double dM_dr = (1.-exp(-y[1]))/2. + r/2. *exp(-y[1])*dy[1];
                                                                                                        return  dM_dr < 1e-18 ; },true);

const integrator::Event FermionBosonStar::Psi_diverging = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return (std::abs(y[3]) > 1.0); }, true);

const integrator::Event FermionBosonStar::phi_negative = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[2] < 0.; });
const integrator::Event FermionBosonStar::phi_positive = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[2] > 0.; });

void FermionBosonStar::set_initial_conditions(const double rho_0, const double phi_0) {
    this->rho_0 = rho_0;
    this->phi_0 = phi_0;
    this->initial_conditions =  vector( {0., 0., phi_0, 0., this->EOS->get_P_from_rho(rho_0, 0.)});
}


int FermionBosonStar::integrate(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts, double r_init, double r_end) {
    return integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  result,  events, intOpts);
}

// find the correct omega-value for a given FBS using bisection in the range [omega_0,omega_1]
// args: FermionBosonStar, vector, min omega, max omega
void FermionBosonStar::bisection(double omega_0, double omega_1, int n_mode, int max_steps, double delta_omega) {

    // values/parameters for bisection
    double omega_mid;
    int n_roots_0, n_roots_1, n_roots_mid;   // number of roots in Phi(r) (number of roots corresponds to the modes of the scalar field)
    //int n_mode = 0;         // number of the mode we want to compute. (mode n has n roots in the scalar field Phi(r))
    int i = 0;

    // variables regarding the integration
    integrator::IntegrationOptions intOpts;
    // define events to check for during integration and put them inside of a std::vector:
    // stop integration if solution diverges:
    std::vector<integrator::Event> events = {phi_negative, phi_positive, Psi_diverging};     // put the events into the event array
    // declare containers to hold the solution of the integration for the upper- (1), lower- (0) and middle (mid) omega
    std::vector<integrator::step> results_0, results_1, results_mid;

    // find initial values for omega min and omega max
    assert(omega_0 < omega_1);  // if the lower omega (omega_0) is larger than the upper omega (omega_1), we stop the code execution

    // set the lower omega and integrate the ODEs:
    this->omega = omega_0;
    int res = this->integrate(results_0, events, intOpts);
    n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    // set the upper omega and integrate the ODEs:
    this->omega = omega_1;
    res = this->integrate(results_1, events, intOpts);
    n_roots_1 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    assert(n_roots_0 != n_roots_1);
    assert(n_roots_0 <= n_mode && n_mode <= n_roots_1);
    //std::cout << "start with omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;

    // find right number of zero crossings (roots) cossesponding to the number of modes (n-th mode => n roots)
    // iterate until the upper and lower omega produce results with one root difference
    while(n_roots_1 - n_roots_0 > 1) {
        omega_mid = (omega_0 + omega_1)/2.;
        //std::cout << "omega_mid = " << omega_mid << " ->";
        this->omega = omega_mid;
        res = this->integrate(results_mid, events, intOpts);
        n_roots_mid = events[0].steps.size() + events[1].steps.size() -1;   // number of roots is number of - to + crossings plus + to - crossings
        //std::cout << " with n_roots = " << n_roots_mid << std::endl;

        if(n_roots_mid == n_roots_0 || n_roots_mid <= n_mode) {
            n_roots_0 = n_roots_mid;
            omega_0 = omega_mid;
            continue;
        }
        if(n_roots_mid == n_roots_1 || n_roots_mid >= n_mode) {
            n_roots_1 = n_roots_mid;
            omega_1 = omega_mid;
            continue;
        }
    }
    //std::cout << "found omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;

    // find right behavior at infty ( Phi(r->infty) = 0 )
    int n_inft_0, n_inft_1, n_inft_mid; // store the sign of Phi at infinity (or at the last r-value)
    this->omega = omega_0; // intOpts.save_intermediate=true;
    res = this->integrate(results_0, events, intOpts);
    n_inft_0 = results_0[results_0.size()-1].second[2] > 0.;    // save if sign(Phi(inf)) is positive or negative

    this->omega = omega_1;
    res = this->integrate(results_1, events, intOpts);
    n_inft_1 = results_1[results_1.size()-1].second[2] > 0.;    // save if sign(Phi(inf)) is positive or negative
    //std::cout << "start with omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << " and omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;

    intOpts.save_intermediate=false;
    while(omega_1 - omega_0 > delta_omega && i < max_steps) { // iterate until accuracy in omega was reached or max number of steps exceeded
        omega_mid = (omega_0 + omega_1)/2.;
        //std::cout << "omega_mid = " << omega_mid << " ->";
        this->omega = omega_mid;
        res = this->integrate(results_mid, events, intOpts);
        n_inft_mid = results_mid[results_mid.size()-1].second[2] > 0.;  // save if sign(Phi(inf)) is positive or negative
        //std::cout << " with n_inft= " << n_inft_mid << std::endl;

        i++;
        // compare the signs of Phi at infinity of the omega-upper, -middle and -lower solution
        // when middle and lower sign are equal, we can move omega_0 to omega_mid
        if(n_inft_mid == n_inft_0) {
            n_inft_0 = n_inft_mid;
            omega_0 = omega_mid;
            continue;
        }
        // when middle and upper sign are equal, we can move omega_1 to omega_mid
        if(n_inft_mid == n_inft_1) {
            n_inft_1 = n_inft_mid;
            omega_1 = omega_mid;
            continue;
        }
    }

    //std::cout << "found omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << " and omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;
    this->omega = omega_0;
}


// uses the bisection method to calculate a FBS solution with fixed rho_c and fixed particle ratio Nb/Nf
// optimized phi_c and omega in the process.
void FermionBosonStar::shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, double omega_0, double omega_1, int n_mode, int max_steps, double delta_omega) {

    // calc the FBS solution once using an initial Phi value
    double my_NbNf;
    double phi_c_init = this->initial_conditions[2];

    while (true) {

        this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
        this->evaluate_model();
        // obtain the ration Nb/Nf
        my_NbNf = this->N_B / this->N_F;
        // check if obtained ratio is above the wanted ratio
        // if yes, we perform a bisection search in the range [0, Phi_initial]
        // if no, we increase Phi_initial by an amount and perform the above steps again

        if (NbNf_ratio < my_NbNf) {
            // the calculated ratio is above the wanted ratio. We can now perform the bisection search!
            break;
        }
        // the wanted ratio is above the calculated ratio. Increase the phi-field for higher ratio Nb/Nf
        if (my_NbNf > 1e-5) {
            phi_c_init = phi_c_init*2.;
        } else{
            phi_c_init = phi_c_init*100.;
        }

        this->set_initial_conditions(this->rho_0, phi_c_init);
        continue;
    }

    // now perform te bisection until the wanted accuracy is reached:
    // define a few local variables:
    double phi_c_0 = 1e-20;
    double phi_c_1 = phi_c_init;
    double phi_c_mid = (phi_c_0 + phi_c_1) / 2.;

    double NbNf_0;
    double NbNf_mid;    // NbNftio of the mid point in phi
    double NbNf_1 = my_NbNf;

    //double my_omega;

    this->set_initial_conditions(this->rho_0, phi_c_0);
    this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
    this->evaluate_model();
    NbNf_0 = this->N_B / this->N_F;

    int i = 0;
    // continue bisection until the wanted accuracy was reached
    while ( (std::abs(NbNf_0 - NbNf_1) > NbNf_accuracy) && (i < max_steps) ) {
        i++;

        phi_c_mid = (phi_c_0 + phi_c_1) / 2.;

        this->set_initial_conditions(this->rho_0, phi_c_mid);
        this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
        //my_omega = this->omega;
        this->evaluate_model();
        // obtain the ration Nb/Nf
        NbNf_mid = this->N_B / this->N_F;

        if (NbNf_mid < NbNf_ratio) {
            // the mid point is below the wanted ratio and we can adjust the lower bound
            NbNf_0 = NbNf_mid;
            phi_c_0 = phi_c_mid;
            continue;
        }
        else if (NbNf_mid > NbNf_ratio) {
            // the mid point is above the wanted ratio and we can adjust the upper bound
            NbNf_1 = NbNf_mid;
            phi_c_1 = phi_c_mid;
            continue;
        }

    }
    // the now obtained omega and phi_c values are now optimized for the wanted Nb/Nf and we can quit the function

}

void FermionBosonStar::evaluate_model() {
    std::vector<integrator::step> results;
    this->evaluate_model(results);
}

void FermionBosonStar::evaluate_model(std::vector<integrator::step>& results, std::string filename) {

    integrator::IntegrationOptions intOpts;
    intOpts.save_intermediate = true;


    std::vector<integrator::Event> events = {FermionBosonStar::M_converged, FermionBosonStar::Psi_diverging};
    //std::vector<integrator::step> results;
    results.clear();

    int res = this->integrate(results, events, intOpts);

    if(!filename.empty()) {
        plotting::save_integration_data(results, {0,1,2,3,4}, {"v", "u", "Phi", "Psi", "P"}, filename);

        #ifdef DEBUG_PLOTTING
        plotting::plot_evolution(results, events, {0,1,2,3,4}, {"v", "u", "Phi", "Psi", "P"}, filename.replace(filename.size()-3, 3, "png"));
        matplotlibcpp::legend(); matplotlibcpp::yscale("log");
        matplotlibcpp::save(filename); matplotlibcpp::close();
        #endif
    }

    // obtain estimate for the total mass:
    // find the minimum in the g_tt component and then compute the total mass M_T at the corresponding index:
    // start iterating through the solution array backwards:
    int min_index_u = results.size()-1;
    auto u_func = [&results](int index) {return results[index].second[1]; };
    double curr_u_min = u_func(results.size()-1);  // last value for the metric component a at r=0

    // find the index of the minimum of the g_rr metric component:
    for (unsigned int i=results.size()-2; i > 0; i--) {
        if (abs(u_func(i)) < curr_u_min) {
            curr_u_min = abs(u_func(i)); // update current minimum
            min_index_u = i;  // update min_index
        }
        else {
            break; // the component is increasing again. that means we have found out minimum
        }
    }

    // find the index of the optimum of the total mass M_T (before it diverges).
    // M_T := r/2 * ( 1 - 1/(a^2) )
    // M_T(r) should i theory be a monotonically increasing function and should therefore have no local minima. Only a global one at r=0
    // Note that this optimum might be at a different position in r than the minimum of the g_rr metric component
    // the optimum can be found at the point where the derivative of M_t with respect to r is minimal:

    int min_index_dMdr = results.size()-1;
    // last value of the total mass M_T:
    auto M_func = [&results](int index) { return results[index].first / 2. * (1. - exp(-results[index].second[1])); };
    double curr_M_last = M_func(results.size()-1),
            curr_M_last_minus1 = M_func(results.size()-2);
    double curr_dMdr_min = (curr_M_last - curr_M_last_minus1) / (results[results.size()-1].first - results[results.size()-2].first);

    for (unsigned int i=results.size()-2; i > 0; i--) {
        curr_M_last = M_func(i);
        curr_M_last_minus1 = M_func(i-1);
        double curr_dMdr_value_i = (curr_M_last - curr_M_last_minus1) / (results[i].first - results[i-1].first);

        if (curr_dMdr_value_i < curr_dMdr_min) {
            curr_dMdr_min = curr_dMdr_value_i; // update current minimum
            min_index_dMdr = i;  // update min_index
            //std::cout << min_index_dMdr << std::endl;
        }
        else {
            break; // the component is increasing again. that means we have found out minimum
        }
    }

    // calculate M_T in where the last local minimum of M_T is:
    double M_T = M_func(min_index_dMdr);

    // std::cout << "min_index_nu: " << min_index_nu << " min_index_M: " << min_index_dMdr << " min_index_phi: " << min_index_phi << " res_size:" << results.size() << std::endl;

    // find the minimum in the phi field before it diverges (to accurately compute the bosonic radius component later):
    // Note: this method will maybe not work well if we consider higher modes of phi!
    int min_index_phi = results.size()-1;
    auto phi_func = [&results](int index) { return results[index].second[2]; };
    double curr_phi_min = abs(phi_func(results.size()-1));  // last value for the phi field a at r=0
    for (unsigned int i=results.size()-2; i >= 0; i--) {
        if (abs(phi_func(i)) < curr_phi_min) {
            curr_phi_min = abs(phi_func(i)); // update current minimum
            min_index_phi = i;  // update min_index
        }
        else {
            break; // the component is increasing again. that means we have found out minimum
        }
    }

    // Extract the results and put them into a usable form to calculate N_B, N_F
    std::vector<double> r(results.size()), N_B_integrand(results.size()), N_F_integrand(results.size());
    vector v;
    double rho, eps;

    for(unsigned int i = 0; i < results.size(); i++) {
        r[i] = results[i].first;
        v = results[i].second;
        N_B_integrand[i] = exp(v[1]/2.) * this->omega *  v[2] * v[2] * r[i] * r[i] * exp(-v[0]/2.);  // get bosonic mass (paricle number) for each r
        this->EOS->callEOS(rho, eps, std::max(0., v[4]));
        N_F_integrand[i] = exp(v[1]/2.) * rho * r[i] * r[i] ;   // get fermionic mass (paricle number) for each r
    }

    // Integrate
    std::vector<double> N_F_integrated, N_B_integrated;
    integrator::cumtrapz(r, N_F_integrand, N_F_integrated);
    integrator::cumtrapz(r, N_B_integrand, N_B_integrated);

    // Find where 99% of N_B,N_F are reached to get the radii
    // we must take the value of the integral *before* the solution diverges! Therefore we cannot just take the last array element
    // but we take the index of the minimum of the metrig g_tt component and the scalar field Phi respectively! This is given by "min_index_*" (see above)
    double N_F = 4.*M_PI* N_F_integrated[min_index_u],
           N_B = 8.*M_PI* N_B_integrated[min_index_phi];

    // first find the index in array where 99% is contained
    // only iterate until the position where the minimum of the metrig g_tt component is (min_index)
    int i_B = 0, i_F = 0;
    unsigned int max_index = std::max(min_index_phi, min_index_u);
    for(unsigned int i = 1; i < max_index; i++) {
        if(N_B_integrated[i] < 0.99*N_B)
            i_B++;
        if(N_F_integrated[i] < 0.99*N_F)
            i_F++;
    }
    // obtain radius from corresponding index
    double R_B = r[i_B], R_F = r[i_F];

    // compute the fermionic radius R_f using the definition where P(R_f)==0:
    // iterate the Pressure-array until we find the first point where the pressure is zero:
    double R_F_0 = 0.0; // fermionic radius where pressure is zero

    // find the first point where the pressure is (approx) zero and take this as the fermionic radius
    for (unsigned i = 1; i < results.size(); i++) {
        if (results[i].second[4] < P_ns_min) {
            R_F_0 = r[i];
            N_F = 4.*M_PI* N_F_integrated[i];
            break;
        }
    }

    //std::cout << "M_T = " << M_T << ", N_B = " << N_B << ", R_B = " << R_B << ", N_F = " << N_F << ", R_F = " << R_F << ", R_F_0 = " << R_F_0 << ", N_B/N_F = " << N_B / N_F << std::endl;
    this->M_T = M_T; this->N_B = N_B; this->N_F = N_F; this->R_B = R_B; this->R_F = R_F; this->R_F_0 = R_F_0;

}

std::ostream& operator<<(std::ostream& os, const FermionBosonStar& fbs) {
    return os   << fbs.M_T                   << " "   // total gravitational mass
                << fbs.rho_0                 << " "   // central density
                << fbs.phi_0                 << " "   // central scalar field
                << fbs.R_F*1.476625061       << " "   // fermionic radius
                << fbs.R_F_0*1.476625061     << " "   // fermionic radius where P(r)=0
                << fbs.N_F                   << " "   // number of fermions
                << fbs.R_B*1.476625061       << " "   // bosonic radius
                << fbs.N_B                   << " "   // number of bosons
                << fbs.N_B / fbs.N_F         << " "   // ratio N_B / N_F
                << fbs.M_T / fbs.R_F_0           << " "   // compactness
                << fbs.omega                 << " "   // omega
                << fbs.mu                    << " "   // mass mu
                << fbs.lambda;      // self-interaction parameter lambda
}


/***********************
 * FermionBosonStarTLN *
 ***********************/

void FermionBosonStarTLN::set_initial_conditions(const double phi_1_0, const double H_0, const double r_init) {
    this->H_0 = H_0;
    this->phi_1_0 = phi_1_0;
    this->initial_conditions =  vector( {0., 0., this->phi_0, 0., this->EOS->get_P_from_rho(this->rho_0, 0.), H_0*r_init*r_init, 2.*H_0*r_init, phi_1_0*pow(r_init,3), 3.*phi_1_0*pow(r_init,2)});
}

vector FermionBosonStarTLN::dy_dt(const double r, const vector& vars) {
    const double v = vars[0], u = vars[1], phi = vars[2], Psi = vars[3];
    double P = vars[4];
    const double H = vars[5],  dH_dr = vars[6],  phi_1 = vars[7], dphi_1_dr = vars[8];

    EquationOfState& myEOS = *(this->EOS);
    const double mu = this->mu; const double lambda = this->lambda; const double omega = this->omega;

    double rho, epsilon, drho_dP, dP_drho;
    if(P <= 0. || P < myEOS.min_P() || P < P_ns_min)  {
        P = 0.; rho = 0.; epsilon = 0., drho_dP = 0.;
    } else {
        myEOS.callEOS(rho, epsilon, P); // change rho and epsilon by reference using EOS member function
        dP_drho = myEOS.dP_drho(rho, epsilon);
        drho_dP = dP_drho > 0. ? 1./dP_drho : 0.;
    }

    vector dy_dr = FermionBosonStar::dy_dt(r, vars); // use equations as given in parent class
    const double dv_dr = dy_dr[0],  du_dr = dy_dr[1], dphi_dr = dy_dr[2], dPsi_dr = dy_dr[3], dP_dr = dy_dr[4];

    const double V = mu*mu*phi*phi + lambda/2.*phi*phi*phi*phi;
    const double dV_deps = mu*mu + lambda*phi*phi;
    const double ddV_deps2 = lambda;

    // additional TLN equations
    const double ddv_dr2 = (1. - exp(u))/r/r
                            + exp(u)*du_dr/r
                            + 8.*M_PI*(omega*omega*exp(u-v)* ( r* (du_dr - dv_dr)*phi*phi + 2.*r*phi*Psi + phi*phi )
                                        + r*exp(u)*du_dr * (P - V) + r*exp(u)* (-2.*phi*Psi*dV_deps +  dP_dr)
                                        + 2.*r*Psi*dPsi_dr + exp(u)*(P - V) + Psi*Psi);

    const double ddH_dr2 = (du_dr/2. - dv_dr/2. - 2./r) * dH_dr
                            + (8.*omega*omega*M_PI*phi*phi*exp(u-v)*(-1.+ drho_dP) + 8.*M_PI *dphi_dr*dphi_dr*(3. + drho_dP)
                                    - ddv_dr2 + du_dr*dv_dr/2. + dv_dr*dv_dr/2. - du_dr/r/2.*(3.+ drho_dP) - dv_dr/r/2.*(7. + drho_dP)
                                    + 6*exp(u)/r/r) * H
                            + (16.*omega*omega*M_PI*phi*exp(u-v)/r*(1. - drho_dP) + 16.*M_PI*phi*exp(u)*dV_deps/r*(1. +drho_dP) - 16.*M_PI* dPsi_dr/r*(3. + drho_dP)
                                    + 8.*M_PI*dphi_dr*du_dr/r *(3. + drho_dP) + 8.*M_PI*dv_dr*dphi_dr/r*(1. - drho_dP) - 32.*M_PI*dphi_dr/r/r*(3. + drho_dP)) * phi_1;

    const double ddphi_1_dr2 = (du_dr - dv_dr)/2.* dphi_1_dr
                                + (omega*omega*r*phi*exp(u-v) - r*dPsi_dr + (r*du_dr/2. + r*dv_dr/2. -2.)*dphi_dr )* H
                                + (-omega*omega*exp(u-v) + 32.*M_PI*Psi*Psi + 2.*phi*phi*exp(u)*ddV_deps2 + exp(u)*dV_deps - du_dr/r/2. + dv_dr/r/2. + 6.*exp(u)/r/r)*phi_1;

    /*std::cout << "r = " << r
                << ", y = " << vars
                << ", dy_dr = " << vector({dy_dr[0], dy_dr[1], dy_dr[2], dy_dr[3], dy_dr[4], dH_dr, ddH_dr2, dphi_1_dr, ddphi_1_dr2})
                << ", ddv_dr2 = " << ddv_dr2
                << ", drho/dP = " << drho_dP
                << std::endl;*/
    return vector({dy_dr[0], dy_dr[1], dy_dr[2], dy_dr[3], dy_dr[4], dH_dr, ddH_dr2, dphi_1_dr, ddphi_1_dr2});
}


void FermionBosonStarTLN::evaluate_model(std::vector<integrator::step>& results, std::string filename) {

    FermionBosonStar::evaluate_model(results, "");


    // add TLN calculation
    // The quantity to compute is y = r H' / H
    // if the fermionic radius larger than the bosonic one, take y = y(R_F_0)
    // if the bosonic radius is larger, find the maxiumum going from the back to the front

    auto y_func = [&results](int index) { return results[index].first * results[index].second[6]/ results[index].second[5]; };
    double y = 0.;

    if(this->R_F_0 > this->R_B) {
        int index_R_F = 0;
        while ( results[index_R_F].first < this->R_F_0 ) index_R_F++;
        // approximate y at R_F_0
        y = y_func(index_R_F-1) +  (y_func(index_R_F) - y_func(index_R_F-1)) / (results[index_R_F].first - results[index_R_F-1].first) * (this->R_F_0 - results[index_R_F-1].first);
        std::cout << "R_F_0 > R_B:  at R_F_0=" << R_F_0 << " y = " << y << std::endl;
    }
    else {
        int max_index_y = results.size()-1;
        // last value of the total mass y:
        double curr_y_last = y_func(results.size()-1);

        for (int i=results.size()-2; i > 0; i--) {
            double curr_y = y_func(i);

            if (curr_y > curr_y_last) {
                curr_y_last = curr_y; // update current maximum
                max_index_y = i;  // update max_index
            }
            else {
                break; // the component is increasing again. that means we have found out minimum
            }
        }
        y = y_func(max_index_y);
        std::cout << "R_B> R_F_0: found max y = " << y << " at " << max_index_y << ", r=" << results[max_index_y].first << std::endl;
    }

    // now that we found y, calculate k2
    double R = std::max(this->R_F_0, this->R_B);
    double C = this->M_T / R; // the compactness

    /* tidal deformability as taken from https://arxiv.org/pdf/0711.2420.pdf */
    double k2 = 8.*pow(C,5)/5. * pow(1.-2.*C, 2)* (2. + 2.*C*(y-1.) - y)
                    / (2.*C*(6. - 3.*y + 3.*C*(5.*y-8.))
                        + 4.*pow(C,3)*(13. - 11.*y + C*(3.*y-2.) + 2.*C*C*(1. + y))
                        + 3.* pow(1. - 2.*C, 2) *(2. - y + 2.*C*(y-1))*log(1.-2.*C));

    /*std::cout << "C = " << C << ", y = " << y << ", k2 = " << k2
                << ", a= " << (2. + 2.*C*(y-1.) - y)

            << std::endl;*/


    this->k2 = k2;

    // add y to results list for easy plotting
    for(int i = 0; i < results.size(); i++) {
        auto s = results[i].second;
        results[i].second = vector({s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], y_func(i)});
    }

    if(!filename.empty()) {
        plotting::save_integration_data(results, {0,1,2,3,4,5,6,7,8,9}, {"a", "alpha", "Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1", "y"}, filename);

        std::vector<integrator::Event> events;
        #ifdef DEBUG_PLOTTING
        plotting::plot_evolution(results, events, {2,3,4,5,6,7,8,9}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1", "y"}, filename.replace(filename.size()-3, 3, "png"), true);
        matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
        matplotlibcpp::save(filename); matplotlibcpp::close();
        #endif
    }
}

std::ostream& operator<<(std::ostream& os, const FermionBosonStarTLN& fbs) {
    return os   << (FermionBosonStar)fbs   << " "   // parent class parameters
                << fbs.k2;  // tidal love number
}

const integrator::Event FermionBosonStarTLN::dphi_1_diverging = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return (std::abs(y[8]) > 1e6); }, true);

const integrator::Event FermionBosonStarTLN::phi_1_negative = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[7] < 0.; });
const integrator::Event FermionBosonStarTLN::phi_1_positive = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[7] > 0.; });


// find the correct phi_1-value for a given FBS using bisection in the range [phi_1_0, phi_1_1]
void FermionBosonStarTLN::bisection_phi_1(double phi_1_0, double phi_1_1, int n_mode, int max_steps, double delta_phi_1) {
    // values/parameters for bisection
    double phi_1_mid;
    int n_roots_0, n_roots_1, n_roots_mid;   // number of roots in phi_1(r) (number of roots corresponds to the modes of the scalar field)
    int i = 0;
    const int index_phi_1 = 7;

    // evaluate model without TLN and calculate values such that we can check consistency at the end
    FermionBosonStar::evaluate_model();

    // variables regarding the integration
    integrator::IntegrationOptions intOpts;
    intOpts.verbose = 1;
    // define events to check for during integration and put them inside of a std::vector:
    // stop integration if solution diverges:
    std::vector<integrator::Event> events = {phi_1_negative, phi_1_positive, dphi_1_diverging};     // put the events into the event array
    // declare containers to hold the solution of the integration for the upper- (1), lower- (0) and middle (mid) phi_1
    std::vector<integrator::step> results_0, results_1, results_mid;

    // find initial values for phi_1 min and phi_1 max
    assert(phi_1_0 < phi_1_1);  // if the lower phi_1 is larger than the upper phi_1

    #ifdef DEBUG_PLOTTING
    intOpts.save_intermediate = true;
    #endif

    // set the lower phi_1 and integrate the ODEs:
    this->set_initial_conditions(phi_1_0, this->H_0);
    int res = this->integrate(results_0, events, intOpts);
    n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    // set the upper phi_1 and integrate the ODEs:
    this->set_initial_conditions(phi_1_1, this->H_0);
    res = this->integrate(results_1, events, intOpts);
    n_roots_1 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    #ifdef DEBUG_PLOTTING
    plotting::plot_evolution(results_0, events, {2,3,4,5,6,7,8}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
    matplotlibcpp::save("test/initial_0.png"); matplotlibcpp::close();
    plotting::plot_evolution(results_1, events, {2,3,4,5,6,7,8}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
    matplotlibcpp::save("test/initial_1.png");
    #endif

    std::cout << "start with phi_1_0 =" << phi_1_0 << " with n_roots=" << n_roots_0 << " and phi_1_1=" << phi_1_1 << " with n_roots=" << n_roots_1 << std::endl;
    assert(n_roots_0 != n_roots_1);
    assert(n_roots_1 <= n_mode && n_mode <= n_roots_0);
    intOpts.save_intermediate = false;

    // find right number of zero crossings (roots) cossesponding to the number of modes (n-th mode => n roots)
    // iterate until the upper and lower phi_1 produce results with one root difference
    while(n_roots_0 - n_roots_1 > 1 && i < max_steps) {
        phi_1_mid = (phi_1_0 + phi_1_1)/2.;
        std::cout << "i=" << i << ": phi_1_mid = " << phi_1_mid << " ->";
        this->set_initial_conditions(phi_1_mid, this->H_0);
        res = this->integrate(results_mid, events, intOpts);
        n_roots_mid = events[0].steps.size() + events[1].steps.size() -1;   // number of roots is number of - to + crossings plus + to - crossings
        std::cout << " with n_roots = " << n_roots_mid << std::endl;
        i++;
        if(n_roots_mid == n_roots_1 || n_roots_mid <= n_mode) {
            n_roots_1 = n_roots_mid;
            phi_1_1 = phi_1_mid;
            continue;
        }
        if(n_roots_mid == n_roots_0 || n_roots_mid >= n_mode) {
            n_roots_0 = n_roots_mid;
            phi_1_0 = phi_1_mid;
            continue;
        }
    }
    std::cout << "after i=" << i << ": found phi_1_0 =" << phi_1_0 << " with n_roots=" << n_roots_0 << " and phi_1_1=" << phi_1_1 << " with n_roots=" << n_roots_1 << std::endl;
    //assert(abs(n_roots_1 - n_roots_0) == 1);

    // find right behavior at infty ( Phi(r->infty) = 0 )
    int n_inft_0, n_inft_1, n_inft_mid; // store the sign of Phi at infinity (or at the last r-value)

    #ifdef DEBUG_PLOTTING
    intOpts.save_intermediate=true;
    #endif
    this->set_initial_conditions(phi_1_0, this->H_0);
    res = this->integrate(results_0, events, intOpts);
    n_inft_0 = results_0[results_0.size()-1].second[index_phi_1] > 0.;    // save if sign(Phi_1(inf)) is positive or negative

    this->set_initial_conditions(phi_1_1, this->H_0);
    res = this->integrate(results_1, events, intOpts);
    n_inft_1 = results_1[results_1.size()-1].second[index_phi_1] > 0.;    // save if sign(Phi_1(inf)) is positive or negative
    std::cout << "start with phi_1_0 =" << phi_1_0 << " with n_inft=" << n_inft_0 << " and phi_1_1=" << phi_1_1 << " with n_inft=" << n_inft_1 << std::endl;

    #ifdef DEBUG_PLOTTING
    plotting::plot_evolution(results_0, events, {2,3,4,5,6,7,8}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
    matplotlibcpp::save("test/intermediate_0.png"); matplotlibcpp::close();
    plotting::plot_evolution(results_1, events, {2,3,4,5,6,7,8}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
    matplotlibcpp::save("test/intermediate_1.png");
    #endif

    intOpts.save_intermediate=false;
    i =0;
    while(phi_1_1 - phi_1_0 > delta_phi_1 && i < max_steps) { // iterate until accuracy in phi_1 was reached or max number of steps exceeded
        phi_1_mid = (phi_1_0 + phi_1_1)/2.;
        std::cout << "i=" << i << ", phi_1_mid = " << phi_1_mid << " ->";
        this->set_initial_conditions(phi_1_mid, this->H_0);
        res = this->integrate(results_mid, events, intOpts);
        n_inft_mid = results_mid[results_mid.size()-1].second[index_phi_1] > 0.;  // save if sign(Phi_1(inf)) is positive or negative
        std::cout << " with n_inft= " << n_inft_mid << std::endl;

        i++;
        // compare the signs of Phi at infinity of the phi_1-upper, -middle and -lower solution
        // when middle and lower sign are equal, we can move phi_1_0 to phi_1_mid
        if(n_inft_mid == n_inft_0) {
            n_inft_0 = n_inft_mid;
            phi_1_0 = phi_1_mid;
            continue;
        }
        // when middle and upper sign are equal, we can move phi_1_1 to phi_1_mid
        if(n_inft_mid == n_inft_1) {
            n_inft_1 = n_inft_mid;
            phi_1_1 = phi_1_mid;
            continue;
        }
    }
    std::cout << "after " << i << " steps found phi_1_0 =" << phi_1_0 << " with n_inft=" << n_inft_0 << " and phi_1_1=" << phi_1_1 << " with n_inft=" << n_inft_1 << std::endl;

    // check for consistency with main equations
    double last_r = results_mid[results_mid.size()-1].first;
    assert(last_r > this->R_F_0  && last_r > this->R_B);

    this->set_initial_conditions(phi_1_0, this->H_0);
}
