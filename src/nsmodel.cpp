#include "nsmodel.hpp"

/* This wrapper function allows to call the class member function dy_dr
 * if the params pointer points to the class instance
 * */
vector NSmodel::dy_dr_static(const double r, const vector& y, const void* params) {
    NSmodel* m = (NSmodel*)params;
    return m->dy_dr(r, y);
}

/* This function calls the integrator for the NSmodel class
 * The dy_dr function is automatically called through the wrapper dy_dr_static class
 * The results are saved in result (only the initial and last step, unless intOpts.save_intermediate=true)
 * A list of events to be tracked during the integration can be passed, but this is optional
 * The initial conditions have to be specified - usually given by NSmodel::get_initial_conditions - but they can be modified
 * The IntegrationOptions will be passed to the integrator
 * The integration starts at r_init and tries to reach r_end
 * The return value is the one given by the integrator, compare integrator::return_reason
 * */
int NSmodel::integrate(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, const vector initial_conditions, integrator::IntegrationOptions intOpts, double r_init, double r_end) const {
    return integrator::RKF45(&(this->dy_dr_static), r_init, initial_conditions, r_end, (void*) this,  result,  events, intOpts);
}

/* This function gives the system of ODEs for the FBS star
 *  for the variables a, alpha, phi, Psi, and P
 *  taken from https://arxiv.org/pdf/2110.11997.pdf
 *
 *  This function is called by the integrator during the integration
 * */
vector FermionBosonStar::dy_dr(const double r, const vector& vars) {

    // rename input & class variables for simpler use
    const double a = vars[0]; const double alpha = vars[1]; const double phi = vars[2]; const double Psi = vars[3]; double P = vars[4];
    EquationOfState& EoS = *(this->EOS);

    // define hydrodynamic quantities
    double rho = 0.;      // restmass density, must be set using EOS
    double epsilon = 0.;  // specific energy denstiy, must be set either through EOS or hydrodynamic relations
    // define potential of the bosonic field
    const double V = mu*mu*phi*phi + lambda/2.*pow(phi, 4);
    const double dV_deps = mu*mu + lambda*phi*phi;

    // apply the EOS
    if(P <= 0. || P < EoS.min_P() || P < P_ns_min)  {
        P = 0.; rho = 0.; epsilon = 0.;
    } else {
        EoS.callEOS(rho, epsilon, P);
    }

    // compute the ODEs:
    double da_dr = 0.5* a *     ( (1.-a*a) / r + 8.*M_PI*r*a*a*( omega*omega*phi*phi/alpha/alpha + V + Psi*Psi/a/a + rho*(1.+epsilon)  ));
    double dalpha_dr = 0.5* alpha * ( (a*a-1.) / r + 8.*M_PI*r*a*a*( omega*omega*phi*phi/alpha/alpha - V  + Psi*Psi/a/a + P ) );
    double dPhi_dr = Psi;
    double dPsi_dr = (-omega*omega*a*a/alpha/alpha + a*a*dV_deps) * phi + (da_dr/a - dalpha_dr/alpha - 2./r) * Psi;
    double dP_dr = -(rho*(1.+epsilon) + P)*dalpha_dr/alpha;

    // write the ODE values into output vector
    return vector({da_dr, dalpha_dr, dPhi_dr, dPsi_dr, dP_dr});
}

/* This event triggers when M is sufficiently converged, i.e. dM_dr < M_T_converged */
const integrator::Event FermionBosonStar::M_converged = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void *params) {
                                                                                                        double dM_dr = ((1. - 1./y[0]/y[0])/2. + r*dy[0]/y[0]/y[0]/y[0]);
                                                                                                        return  dM_dr < M_T_converged ; }, true, "M_converged");

/* This event triggers when Psi starts to diverge, i.e. Psi > 1. */
const integrator::Event FermionBosonStar::Psi_diverging = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (std::abs(y[3]) > 1.0); }, true, "Psi_diverging");

/* This event triggers when Psi becomes positive - why is this here? to force an accurate resolution wherever Psi becomes positive, so that phi_converged might trigger */
const integrator::Event FermionBosonStar::Psi_positive = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (y[3] > 0.0); }, false, "Psi_positive", 1e-8);

/* This event triggers when phi becomes negative, it is used to count the zero crossings of phi */
const integrator::Event FermionBosonStar::phi_negative = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return y[2] < 0.; }, false, "phi_negative");
/* This event triggers when phi becomes positive */
const integrator::Event FermionBosonStar::phi_positive = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return y[2] > 0.; }, false, "phi_positive");

/* This event triggers when phi is sufficiently converged, i.e. phi < PHI_converged and dphi < PHI_converged */
const integrator::Event FermionBosonStar::phi_converged = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (std::abs(y[2])/((FermionBosonStar*)params)->phi_0 < PHI_converged &&
                                                                                                            std::abs(y[3])/((FermionBosonStar*)params)->phi_0/r < PHI_converged); }, true, "phi_converged");
                                                                                                /*{ return (std::abs(y[2])/((FermionBosonStar*)params)->phi_0 < PHI_converged && y[3] > 0. && dy[3] > 0.); }, true, "phi_converged");*/

/* This event triggers when the whole integration (the first five variables - specifically to exclude H, dH in TLN integration) is sufficiently converged */
const integrator::Event FermionBosonStar::integration_converged = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return ublas::norm_inf( dy.sub_range(0, 4))/dr/r < PHI_converged; }, true, "integration_converged");

/* This event triggers when the whole integration (the first five variables - specifically to exclude H, dH in TLN integration) is sufficiently converged */
const integrator::Event FermionBosonStar::P_min_reached = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (y[4] <= std::max(P_ns_min, ((FermionBosonStar*)params)->EOS->min_P())); }, false, "P_min_reached", 1e-5);

/* This function gives the initial conditions for the FBS integration
 * with a_0 = 1,  alpha_0 = 1,  phi_0 = this->phi_0, Psi = 0, and P_0 = EOS(rho_0)
 * at the position r_init (which is assumed to be small enough for this to be valid) */
vector FermionBosonStar::get_initial_conditions(const double r_init) const {
    return vector( {1.0, 1.0, this->phi_0, 0., rho_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(this->rho_0, 0.) : 0.});
}



/* This function integrates the FBS system of equations and stops when phi is sufficiently converged (via the phi_converged event)
 * It then restarts the integration after artifically setting phi = Psi = 0 at that point
 * This allows to integrate the neutron matter and metric components until they converge
 *
 * Only call after omega is set to the corresponding value, otherwise this function is not useful.
 * */
int FermionBosonStar::integrate_and_avoid_phi_divergence(std::vector<integrator::step>& results, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts, std::vector<int> additional_zero_indices, double r_init, double r_end) const {

    results.clear();

    integrator::Event Psi_positive = FermionBosonStar::Psi_positive;
    integrator::Event phi_converged = FermionBosonStar::phi_converged;
    if(this->phi_0 > 0.) // if there is no phi component ignore
        phi_converged.stopping_condition = true;
    events.push_back(Psi_positive);
    events.push_back(phi_converged);

    int res = this->integrate(results, events, this->get_initial_conditions(), intOpts, r_init, r_end);

    if(this->phi_0 <= 0.) // no divergence, so easy peasy return
        return res;

    phi_converged = events[events.size()-1];
    if (!phi_converged.active ) { // the integrator didn't catch the convergence so we have to find it ourselves

        if (!intOpts.save_intermediate) { // can't find it if we don't have the steps
            results.clear();
            intOpts.save_intermediate = true;
            int res = this->integrate(results, events, this->get_initial_conditions(), intOpts, r_init, r_end);
            intOpts.save_intermediate = false;
        }
        // find the minimum in the phi field before it diverges
        int index_phi_converged = 1;

        auto abs_phi_func = [&results] (int index) { return std::abs(results[index].second[2]); };
        std::vector<int> abs_phi_minima({});
        int index_phi_global_min = 0;
        for (unsigned int i=1; i < results.size()-1; i++) {
            if ( abs_phi_func(i) <= abs_phi_func(i-1) && abs_phi_func(i) < abs_phi_func(i+1) )
                abs_phi_minima.push_back(i);
            if (abs_phi_func(i) < abs_phi_func(index_phi_global_min) )
                index_phi_global_min = i;
        }
        index_phi_converged = index_phi_global_min;
        if(abs_phi_minima.size() > 0)
            index_phi_converged = abs_phi_minima[abs_phi_minima.size()-1];

        results.erase(results.begin()+index_phi_converged, results.end()); // remove elements from integration

        for (auto it = events.begin(); it != events.end(); ++it) // set events to active so that they don't trigger again in case they were active at the point of convergence
            it->active = true;
    }

    // now we can restart the integration
    vector initial_conditions = results[results.size()-1].second;
    initial_conditions[2] = 0.; initial_conditions[3] = 0.; // artificially set phi = Psi = 0
    for (auto it = additional_zero_indices.begin(); it != additional_zero_indices.end(); ++it)  // and any other specified
        initial_conditions[*it] = 0.;

    events[events.size()-1].stopping_condition = false;
    std::vector<integrator::step> additional_results;
    events.push_back(FermionBosonStar::integration_converged);
    double last_r = results[results.size()-1].first;
    intOpts.clean_events = false;  // to stop the integrator from clearing the events

    res = this->integrate(additional_results, events, initial_conditions, intOpts, last_r, r_end);

    results.reserve(results.size() + additional_results.size());
    for(unsigned int i = 1; i < additional_results.size(); i++)     // skip the first element to avoid duplicates
        results.push_back(additional_results[i]);

    return res;
}

// find the correct omega-value for a given FBS using bisection in the range [omega_0,omega_1]
// args: FermionBosonStar, vector, min omega, max omega
/* This function tries to find the corresponding omega for a given mu, lambda, rho_0, phi_0
 * such that it is an eigenfrequency with a bisection algorithm.
 * The initial range is given by (omega_0, omega_1). First, the right number of zero crossings for (n_mode) is found
 * with the bisection criterion. Then the algorithm tries to hone in on the value of omega such that phi->0 at infty.
 * This is numerically impossible, as phi always diverges at finite r. This can be used as a bisection criterion, checking
 *  whether phi diverges to positive or negative infty. In this way we converge on an omega such that the integration diverges
 *  as late as possible, within max_steps or until a desired accuracy of delta_omega is reached.
 * If the initial range is insufficient, the function tries to extend it
 * The resulting omega is saved in the class
 * */
int FermionBosonStar::bisection(double omega_0, double omega_1, int n_mode, int max_steps, double delta_omega, int verbose) {

    // values/parameters for bisection
    double omega_mid;
    int n_roots_0, n_roots_1, n_roots_mid;   // number of roots in Phi(r) (number of roots corresponds to the modes of the scalar field)
    int i = 0; // count the steps for comparison with max_step

    // variables for integration
    integrator::IntegrationOptions intOpts;
    std::vector<integrator::Event> events = {phi_negative, phi_positive, Psi_diverging};
    std::vector<integrator::step> results_0, results_1, results_mid;

    if( omega_1 < omega_0)
        std::swap(omega_0, omega_1);

    // if phi_0 = 0 then we don't need a bisection
    if (this->phi_0 == 0.) {
        this->omega = 0.;
        return 0;
    }

    // set the lower omega and integrate the ODEs:
    this->omega = omega_0;
    int res = this->integrate(results_0, events, this->get_initial_conditions(), intOpts);
    n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    // set the upper omega and integrate the ODEs:
    this->omega = omega_1;
    res = this->integrate(results_1, events, this->get_initial_conditions(), intOpts);
    n_roots_1 = events[0].steps.size() + events[1].steps.size() - 1;

    // if the desired number of roots is not given by the initial range adjust the range
    if(n_roots_0 == n_roots_1 || n_roots_0 > n_mode || n_mode > n_roots_1) {
        const int max_tries = 20;
        int tries = 0;
        if (verbose > 0)
            std::cout << "omega range insufficient. adjusting range..." << "\n"
                        << "start with omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;
        // adjust omega_0 if it is too large:
        while (n_roots_0 > n_mode) {
            // set the new lower omega and integrate the ODEs:
            omega_0 *= 0.5;
            this->omega = omega_0;
            if (verbose > 1)
                std::cout << tries << ": omega_0 now= " << this->omega << std::endl;
            int res = this->integrate(results_0, events, this->get_initial_conditions(), intOpts);
            n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1; // number of roots is number of - to + crossings plus + to - crossings
            if(tries > max_tries)
                return -1;
            tries++;
        }
        // adjust omega_1 if it is too small:
        while (n_mode >= n_roots_1) {
            // set the new upper omega and integrate the ODEs:
            omega_1 *= 2.0;
            this->omega = omega_1;
            if (verbose > 1)
                std::cout << tries << ": omega_1 now= " << this->omega << std::endl;
            res = this->integrate(results_1, events, this->get_initial_conditions(), intOpts);
            n_roots_1 = events[0].steps.size() + events[1].steps.size() - 1; // number of roots is number of - to + crossings plus + to - crossings
            if(tries > max_tries)
                return -1;
            tries++;
        }
        if (verbose > 0)
            std::cout << "adjusted omega range successfully with omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;
    }

    if (verbose > 0)
        std::cout << "start with omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;

    // find right number of zero crossings (roots) cossesponding to the number of modes (n-th mode => n roots)
    // iterate until the upper and lower omega produce results with one root difference
    while(n_roots_1 - n_roots_0 > 1) {
        omega_mid = (omega_0 + omega_1)/2.;
        this->omega = omega_mid;
        res = this->integrate(results_mid, events, this->get_initial_conditions(), intOpts);
        n_roots_mid = events[0].steps.size() + events[1].steps.size() -1;   // number of roots is number of - to + crossings plus + to - crossings

        if (verbose> 1)
            std::cout << "omega_mid = " << omega_mid  << " with n_roots = " << n_roots_mid << std::endl;

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
    if (verbose > 0)
        std::cout << "found omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;

    // find right behavior at infty ( Phi(r->infty) = 0 )
    int n_inft_0, n_inft_1, n_inft_mid; // store the sign of Phi at infinity (or at the last r-value)
    this->omega = omega_0;
    res = this->integrate(results_0, events, this->get_initial_conditions(), intOpts);
    n_inft_0 = results_0[results_0.size()-1].second[2] > 0.;    // save if sign(Phi(inf)) is positive or negative

    this->omega = omega_1;
    res = this->integrate(results_1, events, this->get_initial_conditions(), intOpts);
    n_inft_1 = results_1[results_1.size()-1].second[2] > 0.;
    if (verbose > 0)
        std::cout << "start with omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << " and omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;

    // optional test if the event triggers for later integration
    integrator::Event phi_converged = FermionBosonStar::phi_converged;
    phi_converged.stopping_condition = false;
    events.push_back(phi_converged);

    while((omega_1 - omega_0)/omega_0 > delta_omega && i < max_steps) { // iterate until accuracy in omega was reached or max number of steps exceeded
        omega_mid = (omega_0 + omega_1)/2.;
        this->omega = omega_mid;
        res = this->integrate(results_mid, events, this->get_initial_conditions(), intOpts);
        n_inft_mid = results_mid[results_mid.size()-1].second[2] > 0.;  // save if sign(Phi(inf)) is positive or negative

        if (verbose > 1)
            std::cout << "omega_mid = " << omega_mid  << " with n_inft= " << n_inft_mid << std::endl;

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
    // TODO: Check usability
    phi_converged = events[events.size()-1];
    if (phi_converged.steps.size() == 0)
        std::cout << "Warning: phi_converged not hasn't triggered after bisection" << std::endl;

    if (n_inft_1 > 0.)
        this->omega = omega_1;
    else
        this->omega = omega_0;
    if (verbose > 0)
        std::cout << "After " << i << " steps, chose " << omega << " out of omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << ", omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;

    return 0;
}


// uses the bisection method to calculate a FBS solution with fixed rho_c and fixed particle ratio Nb/Nf
// optimized phi_c and omega in the process.
void FermionBosonStar::shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, double omega_0, double omega_1, int n_mode, int max_steps, double delta_omega) {

    // calc the FBS solution once using an initial Phi value
    double my_NbNf;
    double phi_c_init = this->phi_0;

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

        this->phi_0 = phi_c_init;
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

    this->phi_0 = phi_c_0;
    this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
    this->evaluate_model();
    NbNf_0 = this->N_B / this->N_F;

    int i = 0;
    // continue bisection until the wanted accuracy was reached
    while ( (std::abs(NbNf_0 - NbNf_1) > NbNf_accuracy) && (i < max_steps) ) {
        i++;

        phi_c_mid = (phi_c_0 + phi_c_1) / 2.;

        this->phi_0 = phi_c_mid;
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

/* This function takes the result of an integration of the FBS system
 * and calculates the star properties
 * M_T, N_B, N_F, R_B, R_F
 * */
void FermionBosonStar::calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events) {
    const int step_number = results.size();

    bool phi_converged = false;
    if (this->phi_0 > 0. && events.size() > 1 && events[events.size()-2].name == FermionBosonStar::phi_converged.name)
        phi_converged = events[events.size()-2].steps.size() > 0;
    std::cout << "calculate_star_parameters with phi_converged = " << phi_converged << std::endl;

    /* find the minimum in the phi field before it diverges (to accurately compute the bosonic radius component later)
     *  */
    int index_phi_converged = 1;

    if (this->phi_0 > 0.) {
        if (phi_converged) {
            auto phi_converged_event = events[events.size()-2];
            while (results[index_phi_converged].first < phi_converged_event.steps[0].first)
                index_phi_converged++;
        } else { // otherwise implement failsafe
            auto abs_phi_func = [&results] (int index) { return std::abs(results[index].second[2]); };
            std::vector<int> abs_phi_minima({});
            int index_phi_global_min = 0;
            for (int i=1; i < step_number-1; i++) {
                if ( abs_phi_func(i) <= abs_phi_func(i-1) && abs_phi_func(i) < abs_phi_func(i+1) )
                    abs_phi_minima.push_back(i);
                if (abs_phi_func(i) < abs_phi_func(index_phi_global_min) )
                    index_phi_global_min = i;
            }
            index_phi_converged = index_phi_global_min;
            if(abs_phi_minima.size() > 0)
                index_phi_converged = std::max(index_phi_converged, abs_phi_minima[abs_phi_minima.size()-1]);
        }
    }

    /*   M_T
     * find the index of the optimum of the total mass M_T (before it diverges)
     * M_T := r/2 * ( 1 - 1/(a^2) )
     * M_T(r) should i theory be a monotonically increasing function and should therefore have no local minima. Only a global one at r=0
     * Note that this optimum might be at a different position in r than the minimum of the g_rr metric component
     * the optimum can be found at the point where the derivative of M_t with respect to r is minimal*/
    double M_T = 0.;
    auto M_func = [&results](int index) { return results[index].first / 2. * (1. - 1./pow(results[index].second[0], 2)); };
    auto dM_func = [&results, &M_func](int i) { return  (M_func(i+1) - M_func(i))/(results[i+1].first - results[i].first)/ 2.
                                                        + (M_func(i) - M_func(i-1))/(results[i].first - results[i-1].first)/2.; };

    if (phi_converged || this->phi_0 == 0.) {
        M_T = M_func(step_number-1);
    }
    else {
        //std::vector<int> dM_minima;
        int index_dM_global_minimum = results.size()-3;
        for (int i=results.size()-3; i > index_phi_converged; i--) {
            if(dM_func(index_dM_global_minimum) != dM_func(index_dM_global_minimum)) // NaN prevention
                index_dM_global_minimum = i;
            /*if(dM_func(i) < dM_func(i-1) && dM_func(i) < dM_func(i+1)) // search for (true) local minima of dM/dr
                dM_minima.push_back(i);*/
            if(dM_func(i) < dM_func(index_dM_global_minimum)) // and find the global one
                index_dM_global_minimum = i;
        }
        // calculate M_T where the last local minimum of M_T is, if it doesn't exist use the global one:
        int min_index_dMdr;
        /*if(dM_minima.size() > 0) {
            std::cout << "found minimum at " << dM_minima[0] << " with " << dM_func(dM_minima[0]) << std::endl;
            min_index_dMdr = dM_minima[0];	// use the first local minimum in the list as it is the one at the largest radius
            min_index_dMdr = min_index_dMdr < index_dM_global_minimum ? index_dM_global_minimum : min_index_dMdr; // the global minimum is actually to the right of the local one, so it should be better
        }
        else*/
            min_index_dMdr = index_dM_global_minimum;
        std::cout << "chose " << min_index_dMdr << " with global minimum " << index_dM_global_minimum << ": " << dM_func(index_dM_global_minimum) << " with step_num " << step_number << std::endl;

        M_T = M_func(min_index_dMdr);
    }

    // std::cout << "min_index_a: " << min_index_a << " min_index_M: " << min_index_dMdr << " min_index_phi: " << min_index_phi << " res_size:" << results.size() << std::endl;


    /*  N_B, N_F
     *  We need to integrate the particle number densities to obtain N_B, N_F */
    std::vector<double> r(step_number), N_B_integrand(step_number), N_F_integrand(step_number);
    vector v;
    double rho, eps;

    for(unsigned int i = 0; i < results.size(); i++) {
        r[i] = results[i].first;
        v = results[i].second;
        N_B_integrand[i] = 8.*M_PI * v[0] * this->omega *  v[2] * v[2] * r[i] * r[i] / v[1];  // get bosonic mass (paricle number) for each r
        if (v[4] < P_ns_min || v[4] < this->EOS->min_P())
            rho = 0.;
        else
            this->EOS->callEOS(rho, eps, v[4]);
        N_F_integrand[i] = 4.*M_PI * v[0] * rho * r[i] * r[i] ;   // get fermionic mass (paricle number) for each r
    }

    // Integrate
    std::vector<double> N_F_integrated, N_B_integrated;
    integrator::cumtrapz(r, N_F_integrand, N_F_integrated);
    integrator::cumtrapz(r, N_B_integrand, N_B_integrated);

    // Find where 99% of N_B,N_F are reached to get the radii
    double N_F =  N_F_integrated[step_number-1],
           N_B =  N_B_integrated[index_phi_converged /*step_number-1*/];

    // first find the index in array where 99% is contained
    // only iterate until the position where the minimum of the metrig g_tt component is (min_index)
    int i_B = 0, i_F = 0;
    unsigned int max_index = /*std::max(min_index_phi, min_index_a)*/ step_number;
    for(unsigned int i = 1; i < max_index; i++) {
        if(N_B_integrated[i] < 0.99*N_B)
            i_B++;
        if(N_F_integrated[i] < 0.99*N_F)
            i_F++;
    }
    // obtain radius from corresponding index
    double R_B = r[i_B], R_F_0 = r[i_F];

    /*  R_F
     * compute the fermionic radius R_f using the definition where P(R_f)==0
     * iterate the Pressure-array until we find the first point where the pressure is zero */
    double R_F = 0.;
    int index_R_F = 0;
    auto P_func = [&results] (int index) {return results[index].second[4];};

    while( P_func(index_R_F) > std::max(P_ns_min, EOS->min_P())  && index_R_F < step_number-1)
        index_R_F++;
    R_F = r[index_R_F];
    N_F = N_F_integrated[index_R_F];

    //std::cout << "M_T = " << M_T << ", N_B = " << N_B << ", R_B = " << R_B << ", N_F = " << N_F << ", R_F = " << R_F << ", R_F_0 = " << R_F_0 << ", N_B/N_F = " << N_B / N_F << std::endl;
    this->M_T = M_T; this->N_B = N_B; this->N_F = N_F; this->R_B = R_B; this->R_F = R_F; this->R_F_0 = R_F_0;
}

/* Simple wrapper function if the results of the integration are not needed for the user */
void FermionBosonStar::evaluate_model() {
    std::vector<integrator::step> results;
    this->evaluate_model(results);
}

/* This function integrates over the ODE system while avoiding the phi divergence
 *  and then calculates the star properties
 *  Optionally, the output of the integration is saved in the file
 *  Only call if omega is the corresponding eigenfrequency of mu, lambda, rho_0, phi_0 */
void FermionBosonStar::evaluate_model(std::vector<integrator::step>& results, integrator::IntegrationOptions intOpts, std::string filename) {

    intOpts.save_intermediate = true;
    intOpts.verbose = 0;

    std::vector<integrator::Event> events;

    integrator::Event P_min_reached = FermionBosonStar::P_min_reached; // the presence of this event will increase the accuracy around R_F
    if(this->rho_0 > 0.)
        events.push_back(P_min_reached);

    int res = this->integrate_and_avoid_phi_divergence(results, events, intOpts);

    if(!filename.empty()) {
        plotting::save_integration_data(results, {0,1,2,3,4}, {"a", "alpha", "phi", "Psi", "P"}, filename);

        #ifdef DEBUG_PLOTTING
        plotting::plot_evolution(results, events, {0,1,2,3,4}, {"a", "alpha", "phi", "Psi", "P"}, filename.replace(filename.size()-3, 3, "png"));
        matplotlibcpp::legend(); matplotlibcpp::xscale("log"); matplotlibcpp::yscale("log");
        matplotlibcpp::save(filename); matplotlibcpp::close();
        #endif
    }

    this->calculate_star_parameters(results, events);
}

/* Outputs the star parameters and properties */
std::ostream& operator<<(std::ostream& os, const FermionBosonStar& fbs) {
    return os   << fbs.M_T                   << " "   // total gravitational mass
                << fbs.rho_0                 << " "   // central density
                << fbs.phi_0                 << " "   // central scalar field
                << fbs.R_F*1.476625061       << " "   // fermionic radius
                << fbs.N_F                   << " "   // number of fermions
                << fbs.R_B*1.476625061       << " "   // bosonic radius
                << fbs.N_B                   << " "   // number of bosons
                << fbs.N_B / fbs.N_F         << " "   // ratio N_B / N_F
                << fbs.omega                 << " "   // omega
                << fbs.mu                    << " "   // mass mu
                << fbs.lambda;      // self-interaction parameter lambda
}
/* Gives the labels of the values from the output */
std::vector<std::string> FermionBosonStar::labels() {
    return std::vector<std::string> ({"M_T", "rho_0", "phi_0", "R_F", "N_F", "R_B", "N_B", "N_B/N_F", "omega", "mu", "lambda"});
}


/***********************
 * FermionBosonStarTLN *
 ***********************/

/* This function gives the initial conditions for the FBS integration
 * with a_0 = 1,  alpha_0 = 1,  phi_0 = this->phi_0, Psi = 0, P_0 = EOS(rho_0)
 *  H_0 = H_0 * r_0^2,  dH_0 = 2 H_0 r , phi_1_0 = phi_0 r_0^3,  dphi_1_0 = phi_0 3 r_0^2
 * at the position r_0=r_init */
vector FermionBosonStarTLN::get_initial_conditions(const double r_init) const {
    return vector( {1.0, 1.0, this->phi_0, 0., this->rho_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(this->rho_0, 0.) : 0.,
                                                             this->H_0*r_init*r_init, 2.*this->H_0*r_init, this->phi_1_0*pow(r_init,3), 3.*this->phi_1_0*pow(r_init,2)});
}

/* This function gives the system of ODEs for the FBS + TLN star
 *  for the variables a, alpha, phi, Psi, P, H, dH, phi_1, dphi_1
 *
 *  This function is called by the integrator during the integration
 * */
vector FermionBosonStarTLN::dy_dr(const double r, const vector& vars) {
    const double a = vars[0], alpha = vars[1], phi = vars[2], Psi = vars[3];
    double P = vars[4];
    const double H = vars[5],  dH_dr = vars[6],  phi_1 = vars[7], dphi_1_dr = vars[8];

    EquationOfState& myEOS = *(this->EOS);

    double rho, epsilon, drho_dP, dP_drho;
    if(P <= 0. || P < myEOS.min_P() || P < P_ns_min)  {
        P = 0.; rho = 0.; epsilon = 0., drho_dP = 0.;
    } else {
        myEOS.callEOS(rho, epsilon, P);
        dP_drho = rho > myEOS.min_rho() ? myEOS.dP_de(rho, epsilon) : 0.;
        drho_dP = dP_drho > 0. ? 1./dP_drho : 0.;
    }

    vector dy_dr = FermionBosonStar::dy_dr(r, vars); // use equations as given in parent class
    const double da_dr = dy_dr[0],  dalpha_dr = dy_dr[1], dphi_dr = dy_dr[2], dPsi_dr = dy_dr[3], dP_dr = dy_dr[4];

    // The equations for the bosonic field potential
    const double V = mu*mu*phi*phi + lambda/2.*phi*phi*phi*phi;
    const double dV_deps = mu*mu + lambda*phi*phi;
    const double ddV_deps2 = lambda;

    // additional TLN equations
    const double ddalpha_dr2 = 4.*M_PI*omega*omega*(2.*r*phi*phi*a*da_dr + 2.*r*phi*a*a*Psi + phi*phi*a*a)/alpha
                                + ( 4.*M_PI*r*a*a*(-omega*omega*phi*phi/alpha/alpha + P - V + Psi*Psi/a/a)
                                    + (a*a - 1.)/2./r ) * dalpha_dr
                                + ( 4.*M_PI*r* ( 2.*P*a*da_dr - 2.*V*a*da_dr - 2.*phi*a*a*Psi*dV_deps + a*a*dP_dr + 2.*Psi*dPsi_dr)
                                    + 4.*M_PI*a*a*(P - V) + 4.*M_PI*Psi*Psi + a*da_dr/r + (1. - a*a)/2./r/r )*alpha;

    const double ddH_dr2 = (da_dr/a - dalpha_dr/alpha - 2./r) * dH_dr
                            + (8.*omega*omega*M_PI*phi*phi*a*a/alpha/alpha*(-1.+ drho_dP) + 8.*M_PI *dphi_dr*dphi_dr*(3. + drho_dP)
                                    - 2.*ddalpha_dr2/alpha + 2.*dalpha_dr*da_dr/alpha/a + 4.*dalpha_dr*dalpha_dr/alpha/alpha - da_dr/r/a*(3.+ drho_dP) - dalpha_dr/r/alpha*(7. + drho_dP)
                                    + 6*a*a/r/r) * H
                            + (16.*omega*omega*M_PI*phi*a*a/r/alpha/alpha*(1. - drho_dP) + 16.*M_PI*phi*a*a*dV_deps/r*(1. +drho_dP) - 16.*M_PI* dPsi_dr/r*(3. + drho_dP)
                                    + 16.*M_PI*dphi_dr*da_dr/r/a *(3. + drho_dP) + 16.*M_PI*dalpha_dr*dphi_dr/r/alpha*(1. - drho_dP) - 32.*M_PI*dphi_dr/r/r*(3. + drho_dP)) * phi_1;

    const double ddphi_1_dr2 = (da_dr/a - dalpha_dr/alpha)* dphi_1_dr
                                + (omega*omega*r*phi*a*a/alpha/alpha - r*dPsi_dr + (r*da_dr/a + r*dalpha_dr/alpha -2.)*dphi_dr )* H
                                + (-omega*omega*a*a/alpha/alpha + 32.*M_PI*Psi*Psi + 2.*phi*phi*a*a*ddV_deps2 + a*a*dV_deps - da_dr/r/a + dalpha_dr/r/alpha + 6.*a*a/r/r)*phi_1;

    /*std::cout << "r = " << r
                << ", ddalpha_dr2 = " << ddalpha_dr2
                << ", drho/dP = " << drho_dP
                << ", ddH_dr2 = " << ddH_dr2
                << ", dH_dr = " << dH_dr
                << ", H = " << H
                << ". ddphi_1_dr2 = " << ddphi_1_dr2
                << ", dphi_1_dr = " << dphi_1_dr
                << ", phi_1 = " << phi_1
                << std::endl;*/
    return vector({dy_dr[0], dy_dr[1], dy_dr[2], dy_dr[3], dy_dr[4], dH_dr, ddH_dr2, dphi_1_dr, ddphi_1_dr2});
}

/* This function takes the result of an integration of the FBS+TLN system
 * and calculates the star properties
 * lambda_tidal, k2, y_max, R_ext
 * */
void FermionBosonStarTLN::calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events) {

    const int step_number = results.size();
    // calculate parameters for unperturbed star
    //FermionBosonStar::calculate_star_parameters(results, events); // TODO: Check if can be uncommented

    /* The quantity to compute is y = r H' / H
     * if the fermionic radius (much) larger than the bosonic one, take y = y(R_F)
     * if the bosonic radius is larger, find the maxiumum going from the back to the front */
    auto M_func = [&results](int index) { return results[index].first / 2. * (1. - 1./pow(results[index].second[0], 2)); };
    auto y_func = [&results](int index) { return results[index].first * results[index].second[6]/ results[index].second[5]; };
    auto dy_func = [&results, &y_func] (int i) { return (y_func(i+1) - y_func(i))/(results[i+1].first - results[i].first)/2. + (y_func(i) - y_func(i-1))/(results[i].first - results[i-1].first)/2.;  };
    double y = 0., R_ext = 0., M_ext = 0.;

    if(this->R_F > 100.*this->R_B) {
        int index_R_F = 0;
        while ( results[index_R_F].first < this->R_F && index_R_F < step_number-1)
            index_R_F++;
        // approximate y, M_ext at R_F
        y = y_func(index_R_F-1) +  (y_func(index_R_F) - y_func(index_R_F-1)) / (results[index_R_F].first - results[index_R_F-1].first) * (this->R_F - results[index_R_F-1].first);
        M_ext = M_func(index_R_F-1) +  (M_func(index_R_F) - M_func(index_R_F-1)) / (results[index_R_F].first - results[index_R_F-1].first) * (this->R_F - results[index_R_F-1].first);
        R_ext = R_F;
        //std::cout << "R_F > R_B:  at R_F=" << R_F << " y = " << y << std::endl;
    }
    else {
        // to find the starting point see where y actually has a minimum
        int index_bs_radius = 1;
        while(results[index_bs_radius].first < this->R_B/1e3  && (unsigned int)index_bs_radius < results.size()-1)
            index_bs_radius++;
        int index_y_min = index_bs_radius;
        while(y_func(index_y_min) < y_func(index_y_min-1) && (unsigned int)index_y_min < results.size()-1)
            index_y_min++;

        // now look for the local maxima&saddle points of y going from left to right (low r to higher r)
        std::vector<int> indices_maxima;
        int i = index_y_min + 1;
        for( unsigned int i = index_y_min + 1; i < results.size()-2; i++) {
            //std::cout << "i=" << i << ", r= " << results[i].first << ", y = " << y_func(i) << ", dy = " << dy_func(i) << std::endl;
            if(   (y_func(i) > y_func(i -1) && y_func(i) > y_func(i+1) )
                    ||  ( y_func(i) > y_func(i-1) && dy_func(i) < dy_func(i-1) && dy_func(i) < dy_func(i+1)) ) {
                indices_maxima.push_back(i);
            }
            if(y_func(i) < 0.) // here something funky is happening so stop
                break;
        }
        int index_y;
        if(indices_maxima.size() == 0) // if nothing was found just take the last point
            index_y = results.size() - 1;
        else
            index_y = indices_maxima.at(indices_maxima.size()-1);
        y = y_func(index_y);
        R_ext = results[index_y].first;
        M_ext = M_func(index_y); // extract M_ext at the same radius
        //std::cout << "R_B> R_F: found max y = " << y << " at " << index_y << ", R_ext=" << R_ext << std::endl;

    }

    // now that we found y, calculate k2
    double C = M_ext / R_ext; // the compactness at the extraction point

    /* tidal deformability as taken from https://arxiv.org/pdf/0711.2420.pdf */
    double lambda_tidal = 16./15.*pow(this->M_T, 5) * pow(1.-2.*C, 2)* (2. + 2.*C*(y-1.) - y)
                    / (2.*C*(6. - 3.*y + 3.*C*(5.*y-8.))
                        + 4.*pow(C,3)*(13. - 11.*y + C*(3.*y-2.) + 2.*C*C*(1. + y))
                        + 3.* pow(1. - 2.*C, 2) *(2. - y + 2.*C*(y-1))*log(1.-2.*C));
    double k2 = 3./2. * lambda_tidal / pow(R_ext, 5);

    /*std::cout << "C = " << C << ", y = " << y  << ", k2 = " << k2
                << ", a= " << (2. + 2.*C*(y-1.) - y)
            << std::endl;*/

    this->lambda_tidal = lambda_tidal; this->k2 = k2; this->y_max = y; this->R_ext= R_ext;
}

void FermionBosonStarTLN::evaluate_model() {
    std::vector<integrator::step> results;
    this->evaluate_model(results);
}

/* This function integrates over the ODE system while avoiding the phi divergence
 *  and then calculates the star properties
 *  Optionally, the output of the integration is saved in the file
 *  Additionally, y is calculated at every point for debugging purposes */
void FermionBosonStarTLN::evaluate_model(std::vector<integrator::step>& results, std::string filename) {

    integrator::IntegrationOptions intOpts;
    intOpts.save_intermediate = true;

    std::vector<integrator::Event> events;
    results.clear();

    std::vector<int> additional_zero_indices = {7,8};
    int res = this->integrate_and_avoid_phi_divergence(results, events,  intOpts, additional_zero_indices);
    /*std::cout << "M_con " << events[0].active << ", Psi_div " << events[1].active << ", dphi_1_div " << events[2].active << std::endl;
    for(auto it = events[0].steps.begin(); it != events[0].steps.end(); ++it) {
        using namespace integrator;
        std::cout << (integrator::step)*it;
    }*/

    auto y_func = [&results](int index) { return results[index].first * results[index].second[6]/ results[index].second[5]; };
    if(!filename.empty()) {
        // add y to results list for easy plotting
        for(unsigned int i = 0; i < results.size(); i++) {
            auto s = results[i].second;
            results[i].second = vector({s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], y_func(i)});
        }
        plotting::save_integration_data(results, {0,1,2,3,4,5,6,7,8,9}, {"a", "alpha", "phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1", "y"}, filename);

        std::vector<integrator::Event> events;
        #ifdef DEBUG_PLOTTING
        plotting::plot_evolution(results, events, {2,3,4,5,6,7,8,9}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1", "y"}, filename.replace(filename.size()-3, 3, "png"));
        matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
        matplotlibcpp::save(filename); matplotlibcpp::close();
        #endif
    }

    this->calculate_star_parameters(results, events);
}

/* This calls the parent class output and adds additional values */
std::ostream& operator<<(std::ostream& os, const FermionBosonStarTLN& fbs) {
    return os   << (FermionBosonStar)fbs   << " "   // parent class parameters
                << fbs.k2                  << " "   // tidal love number
                << fbs.lambda_tidal        << " "   // dimensionfull tidal love number
                << fbs.phi_1_0             << " "   // phi_1 initial value
                << fbs.H_0                          // H initial value
                ;
}

std::vector<std::string> FermionBosonStarTLN::labels() {
    auto l = FermionBosonStar::labels();
    l.push_back("k2"); l.push_back("lambda_tidal"); l.push_back("phi_1_0"); l.push_back("H_0");
    return l;
}

/* This event triggers when dphi_1 is diverging i.e. dphi_1 > 1e6 */
const integrator::Event FermionBosonStarTLN::dphi_1_diverging = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return (std::abs(y[8]) > 1e6); }, true);

/* This event triggers when phi_1 is negative, used for zero counting */
const integrator::Event FermionBosonStarTLN::phi_1_negative = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[7] < 0.; });
/* This event triggers when phi_1 is positive  */
const integrator::Event FermionBosonStarTLN::phi_1_positive = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[7] > 0.; });


/* This function finds the corresponding phi_1_0 inside the range (phi_1_0_l, phi_1_0_r)
 *  such that phi_1 adheres to the boundary conditions phi_1->0 at infty
 *  via a bisection algorithm, similarly to omega
 * This algorithm does not adjust the range
 * The result depends on H_0, do not change afterwards */
int FermionBosonStarTLN::bisection_phi_1(double phi_1_0_l, double phi_1_0_r, int n_mode, int max_steps, double delta_phi_1, int verbose) {

    double phi_1_0_mid;
    int n_roots_0, n_roots_1, n_roots_mid;   // number of roots in phi_1(r) (number of roots corresponds to the modes of the scalar field)
    int i = 0;
    const int index_phi_1 = 7;

    // variables regarding the integration
    integrator::IntegrationOptions intOpts;
    intOpts.verbose = 0;
    std::vector<integrator::Event> events = {phi_1_negative, phi_1_positive, dphi_1_diverging};
    std::vector<integrator::step> results_0, results_1, results_mid;

    if(phi_1_0_r < phi_1_0_l)
        std::swap(phi_1_0_l, phi_1_0_r);

    // if phi_0 = 0 then we don't need a bisection
    if (this->phi_0 == 0.) {
        this->phi_1_0 = 0.;
        return 0;
    }

    #ifdef DEBUG_PLOTTING
    intOpts.save_intermediate = true;
    #endif

    // set the lower phi_1 and integrate the ODEs:
    this->phi_1_0 = phi_1_0_l;
    int res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_0, events, intOpts);
    n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    // set the upper phi_1 and integrate the ODEs:
    this->phi_1_0 = phi_1_0_r;
    res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_1, events, intOpts);
    n_roots_1 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    #ifdef DEBUG_PLOTTING
    plotting::plot_evolution(results_0, events, {2,3,4,5,6,7,8}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
    matplotlibcpp::save("test/initial_0.png"); matplotlibcpp::close();
    plotting::plot_evolution(results_1, events, {2,3,4,5,6,7,8}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
    matplotlibcpp::save("test/initial_1.png");
    #endif

    if(n_roots_0 == n_roots_1 || n_roots_1 > n_mode || n_mode > n_roots_0)
        return -1;

    if (verbose > 0)
        std::cout << "start with phi_1_0_l =" << phi_1_0_l << " with n_roots=" << n_roots_0 << " and phi_1_r=" << phi_1_0_r << " with n_roots=" << n_roots_1 << std::endl;
    intOpts.save_intermediate = false;

    /* find right number of zero crossings (roots) cossesponding to the number of modes (n-th mode => n roots)
     * iterate until the upper and lower phi_1 produce results with one root difference */
    while(n_roots_0 - n_roots_1 > 1 && i < max_steps) {
        phi_1_0_mid = (phi_1_0_l + phi_1_0_r)/2.;
        this->phi_1_0 = phi_1_0_mid;
        res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_mid, events, intOpts);
        n_roots_mid = events[0].steps.size() + events[1].steps.size() -1;   // number of roots is number of - to + crossings plus + to - crossings
        if (verbose > 1)
            std::cout << "i=" << i << ": phi_1_0_mid = " << phi_1_0_mid << " with n_roots = " << n_roots_mid << std::endl;
        i++;
        if(n_roots_mid == n_roots_1 || n_roots_mid <= n_mode) {
            n_roots_1 = n_roots_mid;
            phi_1_0_r = phi_1_0_mid;
            continue;
        }
        if(n_roots_mid == n_roots_0 || n_roots_mid >= n_mode) {
            n_roots_0 = n_roots_mid;
            phi_1_0_l = phi_1_0_mid;
            continue;
        }
    }
    if (verbose > 0)
        std::cout << "after i=" << i << ": found phi_1_0_l =" << phi_1_0_l << " with n_roots=" << n_roots_0 << " and phi_1_0_r =" << phi_1_0_r << " with n_roots=" << n_roots_1 << std::endl;
    if(abs(n_roots_1 - n_roots_0) != 1) // number of roots does no match, we can't continue
        return -1;

    // find right behavior at infty ( Phi(r->infty) = 0 )
    int n_inft_0, n_inft_1, n_inft_mid;

    #ifdef DEBUG_PLOTTING
    intOpts.save_intermediate=true;
    #endif
    this->phi_1_0 = phi_1_0_l;
    res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_0, events, intOpts);
    n_inft_0 = results_0[results_0.size()-1].second[index_phi_1] > 0.;    // save if sign(Phi_1(inf)) is positive or negative

    this->phi_1_0 = phi_1_0_r;
    res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_1, events, intOpts);
    n_inft_1 = results_1[results_1.size()-1].second[index_phi_1] > 0.;

    if (verbose > 0)
        std::cout << "start with phi_1_0_l =" << phi_1_0_l << " with n_inft=" << n_inft_0 << " and phi_1_0_r =" << phi_1_0_r << " with n_inft=" << n_inft_1 << std::endl;

    #ifdef DEBUG_PLOTTING
    plotting::plot_evolution(results_0, events, {2,3,4,5,6,7,8}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
    matplotlibcpp::save("test/intermediate_0.png"); matplotlibcpp::close();
    plotting::plot_evolution(results_1, events, {2,3,4,5,6,7,8}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
    matplotlibcpp::save("test/intermediate_1.png");matplotlibcpp::close();
    #endif

    intOpts.save_intermediate=false;
    i =0;
    /* iterate until accuracy in phi_1 was reached or max number of steps exceeded */

    while( (phi_1_0_r - phi_1_0_l)/phi_1_0_l > delta_phi_1 && i < max_steps) {
        phi_1_0_mid = (phi_1_0_l + phi_1_0_r)/2.;
        this->phi_1_0 = phi_1_0_mid;
        res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_mid, events, intOpts);
        n_inft_mid = results_mid[results_mid.size()-1].second[index_phi_1] > 0.;  // save if sign(Phi_1(inf)) is positive or negative
        if (verbose > 1)
            std::cout << "i=" << i << ", phi_1_0_mid = " << phi_1_0_mid << " with n_inft= " << n_inft_mid << std::endl;

        i++;
        // compare the signs of Phi at infinity of the phi_1-upper, -middle and -lower solution
        // when middle and lower sign are equal, we can move phi_1_0_l to phi_1_0_mid
        if(n_inft_mid == n_inft_0) {
            n_inft_0 = n_inft_mid;
            phi_1_0_l = phi_1_0_mid;
            continue;
        }
        // when middle and upper sign are equal, we can move phi_1_0_r to phi_1_0_mid
        if(n_inft_mid == n_inft_1) {
            n_inft_1 = n_inft_mid;
            phi_1_0_r = phi_1_0_mid;
            continue;
        }
    }

    #ifdef DEBUG_PLOTTING
    intOpts.save_intermediate=true;
    this->phi_1_0 = phi_1_0_l;
    res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_0, events, intOpts);

    plotting::plot_evolution(results_0, events, {2,3,4,5,6,7,8}, {"Phi", "Psi", "P", "H", "dH", "phi_1", "dphi_1"});
    matplotlibcpp::legend(); matplotlibcpp::yscale("log"); matplotlibcpp::xscale("log");
    matplotlibcpp::save("test/final.png"); matplotlibcpp::close();
    #endif

    double last_r = results_mid[results_mid.size()-1].first;
    if (verbose > 0)
        std::cout << "after " << i << " steps found phi_1_0_l =" << phi_1_0_l << " with n_inft=" << n_inft_0 << " and phi_1_0_r=" << phi_1_0_r << " with n_inft=" << n_inft_1
                    << "\n  last_r = " << last_r << " vs R_F_0 = " << this->R_F_0 << "and R_B = " << this->R_B << std::endl;
    //assert(last_r > this->R_F_0  && last_r > this->R_B);

    this->phi_1_0 = phi_1_0_l;
    return 0;
}

/* This function integrates the FBS+TLN system of equations and stops when phi is sufficiently converged (via the phi_converged event)
 * It then restarts the integration after artifically setting phi = Psi = phi_1 = dphi_1 = 0 at that point
 * This allows to integrate the neutron matter and metric components until they converge
 * */
/*
int FermionBosonStarTLN::integrate_and_avoid_phi_divergence(std::vector<integrator::step>& results, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts, double r_init, double r_end) const {

    results.clear();

    integrator::Event Psi_positive = FermionBosonStar::Psi_positive;
    integrator::Event phi_converged = FermionBosonStar::phi_converged;
    if(this->phi_0 > 0.) // if there is no phi component ignore
        phi_converged.stopping_condition = true;
    events.push_back(Psi_positive);
    events.push_back(phi_converged);

    int res = this->integrate(results, events, this->get_initial_conditions(), intOpts, r_init, r_end);

    phi_converged = events[events.size()-1];
    if(res == integrator::event_stopping_condition && phi_converged.active) {
        //std::cout << "phi has converged before P=0 -> restart integration" << std::endl;
        vector initial_conditions = results[results.size()-1].second;
        initial_conditions[2] = 0.; initial_conditions[3] = 0.; // artificially set phi = Psi = 0
        initial_conditions[7] = 0.; initial_conditions[8] = 0.; // set phi_1 = dphi_1 = 0

        events[events.size()-1].stopping_condition = false;
        std::vector<integrator::step> additional_results;
        events.push_back(FermionBosonStar::integration_converged);
        double last_r = results[results.size()-1].first;
        intOpts.clean_events = false;  // to stop the integrator from clearing the events

        res = this->integrate(additional_results, events, initial_conditions, intOpts, last_r, r_end);

        results.reserve(results.size() + additional_results.size());
        for(unsigned int i = 1; i < additional_results.size(); i++)     // skip the first element to avoid duplicates
            results.push_back(additional_results[i]);
    }
    return res;
}*/
