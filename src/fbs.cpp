#include "fbs.hpp"

/*   Events   */
/* This event triggers when M is sufficiently converged, i.e. dM_dr < M_T_converged
 *  (not in use anymore) */
const integrator::Event FermionBosonStar::M_converged = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void *params) {
                                                                                                        const NUMERIC dM_dr = ((1._num - 1._num/y[0]/y[0])/2._num + r*dy[0]/y[0]/y[0]/y[0]);
                                                                                                        return  dM_dr < M_T_converged ; }, true, "M_converged");

/* This event triggers when Psi starts to diverge, i.e. Psi > 1. */
const integrator::Event FermionBosonStar::Psi_diverging = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (std::abs(y[3]) > 1._num); }, true, "Psi_diverging");

/* This event triggers when Psi becomes positive - why is this here? to force an accurate resolution wherever Psi becomes positive, so that phi_converged might trigger */
const integrator::Event FermionBosonStar::Psi_positive = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (y[3] > 0._num); }, false, "Psi_positive", 1e-8_num);

/* This event triggers when phi becomes negative, it is used to count the zero crossings of phi */
const integrator::Event FermionBosonStar::phi_negative = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return y[2] < 0._num; }, false, "phi_negative");
/* This event triggers when phi becomes positive */
const integrator::Event FermionBosonStar::phi_positive = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return y[2] > 0._num; }, false, "phi_positive");

/* This event triggers when phi is sufficiently converged, i.e. phi < PHI_converged and dphi < PHI_converged */
const integrator::Event FermionBosonStar::phi_converged = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (std::abs(y[2])/((FermionBosonStar*)params)->phi_0 < PHI_converged &&
                                                                                                            std::abs(y[3])/((FermionBosonStar*)params)->phi_0/r < PHI_converged); }, true, "phi_converged");
                                                                                                /*{ return (std::abs(y[2])/((FermionBosonStar*)params)->phi_0 < PHI_converged && y[3] > 0. && dy[3] > 0.); }, true, "phi_converged");*/

/* This event triggers when the whole integration (the first five variables - specifically to exclude H, dH in TLN integration) is sufficiently converged */
const integrator::Event FermionBosonStar::integration_converged = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return ublas::norm_inf( dy.sub_range(0, 4))/dr/r < INT_converged; }, true, "integration_converged");

/* This event triggers when the neutron has has converged  */
const integrator::Event FermionBosonStar::P_min_reached = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (y[4] <= std::max(P_ns_min, ((FermionBosonStar*)params)->EOS->min_P())); }, false, "P_min_reached", 1e-5_num);


/* This function gives the system of ODEs for the FBS star
 *  for the variables a, alpha, phi, Psi, and P
 *  taken from https://arxiv.org/pdf/2110.11997.pdf
 *
 *  This function is called by the integrator during the integration
 * */
vector FermionBosonStar::dy_dr(const NUMERIC r, const vector& vars) const {

    // rename input & class variables for simpler use
    const NUMERIC a = vars[0]; const NUMERIC alpha = vars[1]; const NUMERIC phi = vars[2]; const NUMERIC Psi = vars[3]; NUMERIC P = vars[4];
    EquationOfState& EoS = *(this->EOS);

    // define hydrodynamic quantities
    NUMERIC etot = 0.;	// total energy density. Must be computes from P using the EoS
    // define potential of the bosonic field
    const NUMERIC V = mu*mu*phi*phi + lambda/2.*pow(phi, 4);
    const NUMERIC dV_deps = mu*mu + lambda*phi*phi;

    // apply the EOS
    if(P <= 0._num || P < EoS.min_P() || P < P_ns_min)  {
        P = 0._num; etot = 0._num;
    } else {
        etot = EoS.get_e_from_P(P);
    }

    // compute the ODEs:
    NUMERIC da_dr = 0.5_num* a *     ( (1._num-a*a) / r + 8._num*M_PI*r*a*a*( omega*omega*phi*phi/alpha/alpha + V + Psi*Psi/a/a + etot  ));
    NUMERIC dalpha_dr = 0.5_num* alpha * ( (a*a-1._num) / r + 8._num*M_PI*r*a*a*( omega*omega*phi*phi/alpha/alpha - V  + Psi*Psi/a/a + P ) );
    NUMERIC dPhi_dr = Psi;
    NUMERIC dPsi_dr = (-omega*omega*a*a/alpha/alpha + a*a*dV_deps) * phi + (da_dr/a - dalpha_dr/alpha - 2._num/r) * Psi;
    NUMERIC dP_dr = -(etot + P)*dalpha_dr/alpha;

    // write the ODE values into output vector
    return vector({da_dr, dalpha_dr, dPhi_dr, dPsi_dr, dP_dr});
}


/* This function gives the initial conditions for the FBS integration
 * with a_0 = 1,  alpha_0 = 1,  phi_0 = this->phi_0, Psi = 0, and P_0 = EOS(rho_0)
 * at the position r_init (which is assumed to be small enough for this to be valid) */
vector FermionBosonStar::get_initial_conditions(NUMERIC r_init) const {
    return vector( {1._num, 1._num, this->phi_0, 0._num, rho_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(this->rho_0, 0._num) : 0._num});
}


/* This function performs an integration and tries to find R_B_0 while integrating, where the bosonic component has sufficiently converged
 *  and we can set it to zero. The convergence criterion is given by the phi_converged event. The function returns R_B_0 by reference.
 * If the condition is not fulfilled, the "convergence" can be forced (with the force boolean) by finding the last minimum of the phi field.
 * */
int FermionBosonStar::find_bosonic_convergence(std::vector<integrator::step>& results, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts, NUMERIC& R_B_0, bool force, NUMERIC r_init, NUMERIC r_end) const {

    if(this->phi_0 <= 0._num)
        return -1;

    integrator::Event Psi_positive = FermionBosonStar::Psi_positive;
    integrator::Event phi_converged = FermionBosonStar::phi_converged;
    phi_converged.stopping_condition = true;
    events.push_back(Psi_positive); // the presence of this event increases the accuracy around zero crossings
    events.push_back(phi_converged);

    int res = this->integrate(results, events, this->get_initial_conditions(), intOpts, r_init, r_end);

    phi_converged = events[events.size()-1];
    events.pop_back(); // Remove events that we added
    events.pop_back();

    if (res == integrator::event_stopping_condition && !phi_converged.active) // in this case another event triggered the stop so just return
        return res;

    if (!phi_converged.active ) { // the integrator didn't catch the convergence so we have to find it ourselves

        if (!intOpts.save_intermediate) { // can't find it if we don't have the steps
            results.clear();
            intOpts.save_intermediate = true;
            res = this->integrate(results, events, this->get_initial_conditions(), intOpts, r_init, r_end);
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

        // maybe the event didn't trigger because psi was too large?
        if(!force) {
            auto y_at_phi_converged = results[index_phi_converged].second;
            vector dy_at_phi_converged = this->dy_dr(results[index_phi_converged].first, y_at_phi_converged);
            y_at_phi_converged[3] =0._num;
            if (!phi_converged.condition(results[index_phi_converged].first, 0._num, y_at_phi_converged, dy_at_phi_converged, (const void*)this)) {
                return res; // no, phi doesn't get close to zero so we shouldn't artifically set it so
            }
        }
        // found a convergence for phi so restart at that point
        results.erase(results.begin()+index_phi_converged, results.end()); // remove elements from integration

        for (auto it = events.begin(); it != events.end(); ++it) // set events to active so that they don't trigger again in case they were active at the point of convergence
            it->active = true;
    }

    R_B_0 = results[results.size()-1].first;
    return res;
}

/* This function integrates the FBS system of equations and stops when phi is sufficiently converged (via the phi_converged event)
 * This convergence is first found with the find_bosonic_convergence function, then R_B_0 is set so that subsequent calls don't have
 *  to find it again.
 * It then restarts the integration after artifically setting phi = Psi = 0 at that point. Additional values
 *  can be set to zero with the additional_zero_indices vector, i.e. phi_1 and dphi_1.
 * This allows to integrate the neutron matter and metric components until they converge as well.
 * If the convergence criterion cannot be fulfilled, it can be forced with the force flag, see find_bosonic_convergence
 *
 * Only call after omega is set to the corresponding value, otherwise this function is not useful.
 * */
int FermionBosonStar::integrate_and_avoid_phi_divergence(std::vector<integrator::step>& results, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts, bool force, std::vector<int> additional_zero_indices, NUMERIC r_init, NUMERIC r_end)  {

    results.clear();
    int res;

    if(this->phi_0 <= 0._num) { // no divergence, so easy peasy return
        res = this->integrate(results, events, this->get_initial_conditions(), intOpts, r_init, r_end);
        return res;
    }

    // integrate to R_B_0
    if(this->R_B_0 == 0._num) { // find it first
        res = this->find_bosonic_convergence(results, events, intOpts, this->R_B_0, force, r_init, r_end);
        if (this->R_B_0 == 0._num) // the algorithm returned without finding R_B_0 so stop integration here (either divergence or other event triggered)
            return res;
        // std::cout << "found R_B_0=" << R_B_0 << std::endl;
    }
    else  {
        res = this->integrate(results, events, this->get_initial_conditions(), intOpts, r_init, this->R_B_0);
        if(res != integrator::endpoint_reached) // e.g. if some other event triggered
            return res;
    }

    // now we can restart the integration
    vector initial_conditions = results[results.size()-1].second;
    initial_conditions[2] = 0._num; initial_conditions[3] = 0._num; // artificially set phi = Psi = 0
    for (auto it = additional_zero_indices.begin(); it != additional_zero_indices.end(); ++it)  // and any other specified, i.e. phi_1, dphi_1
        initial_conditions[*it] = 0._num;

    std::vector<integrator::step> additional_results;
    events.push_back(FermionBosonStar::integration_converged);
    NUMERIC last_r = results[results.size()-1].first;
    intOpts.clean_events = false;  // to stop the integrator from clearing the events

    res = this->integrate(additional_results, events, initial_conditions, intOpts, last_r, r_end);

    // add the results together into one vector
    results.reserve(results.size() + additional_results.size());
    for(unsigned int i = 1; i < additional_results.size(); i++)     // skip the first element to avoid duplicates
        results.push_back(additional_results[i]);

    events.pop_back();
    return res;
}

int FermionBosonStar::bisection_expand_range(NUMERIC& omega_0, NUMERIC& omega_1, int n_mode, int& n_roots_0, int& n_roots_1, int verbose) {

    integrator::IntegrationOptions intOpts;
    intOpts.verbose = verbose - 1;
    std::vector<integrator::Event> events = {phi_negative, phi_positive, Psi_diverging};
    std::vector<integrator::step> results_0, results_1;

    // if the desired number of roots is not given by the initial range adjust the range
    const int max_tries = 20;
    int tries = 0;
    if (verbose > 0)
        std::cout << "omega range insufficient. adjusting range..." << "\n"
                    << "start with omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;

    while (n_roots_0 > n_mode) {
        // set the new lower omega and integrate the ODEs:
        omega_0 *= 0.5_num;
        this->omega = omega_0;
        if (verbose > 1)
            std::cout << tries << ": omega_0 now= " << this->omega << std::endl;
        this->integrate(results_0, events, this->get_initial_conditions(), intOpts);
        n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1; // number of roots is number of - to + crossings plus + to - crossings
        if(tries > max_tries)
            return -1;
        tries++;
    }

    while (n_mode >= n_roots_1) {
        // set the new upper omega and integrate the ODEs:
        omega_1 *= 2._num;
        this->omega = omega_1;
        if (verbose > 1)
            std::cout << tries << ": omega_1 now= " << this->omega << std::endl;
        this->integrate(results_1, events, this->get_initial_conditions(), intOpts);
        n_roots_1 = events[0].steps.size() + events[1].steps.size() - 1; // number of roots is number of - to + crossings plus + to - crossings
        if(tries > max_tries)
            return -1;
        tries++;
    }
    if (verbose > 0)
        std::cout << "adjusted omega range successfully with omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;
    return 0;
}


int FermionBosonStar::bisection_find_mode(NUMERIC& omega_0, NUMERIC& omega_1, int n_mode, int max_steps, int verbose) {

    NUMERIC omega_mid;
    int n_roots_0, n_roots_1, n_roots_mid;   // number of roots in Phi(r) (number of roots corresponds to the modes of the scalar field)
    int res;
    int steps = 0;

    integrator::IntegrationOptions intOpts;
    intOpts.verbose = verbose - 1;
    std::vector<integrator::Event> events = {phi_negative, phi_positive, Psi_diverging};
    std::vector<integrator::step> results_0, results_1, results_mid;

    // set the lower omega and integrate the ODEs:
    this->omega = omega_0;
    res = this->integrate(results_0, events, this->get_initial_conditions(), intOpts);
    n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    // set the upper omega and integrate the ODEs:
    this->omega = omega_1;
    res = this->integrate(results_1, events, this->get_initial_conditions(), intOpts);
    n_roots_1 = events[0].steps.size() + events[1].steps.size() - 1;

    if(n_roots_0 == n_roots_1 || n_roots_0 > n_mode || n_mode > n_roots_1) {
        res = this->bisection_expand_range(omega_0, omega_1, n_mode, n_roots_0, n_roots_1, verbose);
        if(res)
            return res;
    }

    if (verbose > 0)
        std::cout << "start with omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;

    // find right number of zero crossings (roots) cossesponding to the number of modes (n-th mode => n roots)
    // iterate until the upper and lower omega produce results with one root difference
    while(n_roots_1 - n_roots_0 > 1 && steps < max_steps) {
        omega_mid = (omega_0 + omega_1)/2._num;
        if (omega_mid == this->omega) { // if this happens we reached floting point accuracy limits
            if (verbose > -1) {
                std::cout << "floating point accuracy reached while searching root at iteration step =" << steps << std::endl;
				std::cout << "no suitable pair of omegas found  omega_0 = " << omega_0 << " with n_roots_0 = " << n_roots_0 << " and omega_1 = " << omega_1 << " with n_roots_1 = " << n_roots_1 << std::endl;
        	}
            break;
        }
        this->omega = omega_mid;
        res = this->integrate(results_mid, events, this->get_initial_conditions(), intOpts);
        n_roots_mid = events[0].steps.size() + events[1].steps.size() -1;   // number of roots is number of - to + crossings plus + to - crossings

        if (verbose> 1)
            std::cout << steps << ": omega_mid = " << omega_mid  << " with n_roots = " << n_roots_mid << std::endl;
		
		// failsave-mechanism if the wanted root is skipped (the bisection may skip multiple roots at once). 
		// In some circumstances it even will skip all roots and then the wanted root is outside the range of the bisection!
		// if(n_roots_mid == n_roots_0 || n_roots_mid == n_roots_1) {
		// 	std::cout << steps << ": omega_mid = " << omega_mid  << " with n_roots = " << n_roots_mid << std::endl;
		// 	std::cout << "shit" << std::endl;
		// }

		// bisection: update the upper and lower bounds:
        if(n_roots_mid == n_roots_0 || n_roots_mid <= n_mode) {
            n_roots_0 = n_roots_mid;
            omega_0 = omega_mid;
        }
        else if(n_roots_mid == n_roots_1 || n_roots_mid >= n_mode) {
            n_roots_1 = n_roots_mid;
            omega_1 = omega_mid;
        }
		// if (verbose > -1) {
		// 		std::cout << "at iteration step = " << steps << ", omega_0 = " << omega_0 << " with n_roots_0 = " << n_roots_0 << " and omega_1 = " << omega_1 << " with n_roots_1 = " << n_roots_1 << std::endl;
        // 	}
        steps++;
    }
    if(abs(n_roots_1 - n_roots_0) != 1) { // number of roots does no match, we can't continue
        if (verbose > 0)
            std::cout << "no suitable pair of omegas found  omega_0 = " << omega_0 << " with n_roots_0 = " << n_roots_0 << " and omega_1 = " << omega_1 << " with n_roots_1 = " << n_roots_1 << std::endl;
        return -1;
    }
    if (verbose > 0)
        std::cout << "found omega_0 =" << omega_0 << " with n_roots=" << n_roots_0 << " and omega_1=" << omega_1 << " with n_roots=" << n_roots_1 << std::endl;

    return 0.;
}

int FermionBosonStar::bisection_converge_through_infty_behavior(NUMERIC omega_0, NUMERIC omega_1, int n_mode, int max_steps, NUMERIC delta_omega, int verbose) {

    NUMERIC omega_mid;
    int res, res_1;
    int steps = 0;

    integrator::IntegrationOptions intOpts;
    intOpts.verbose = verbose - 1;
    std::vector<integrator::Event> events = {Psi_diverging};
    std::vector<integrator::step> results_0, results_1, results_mid;

    // find right behavior at infty ( Phi(r->infty) = 0 )
    int n_inft_0, n_inft_1, n_inft_mid; // store the sign of Phi at infinity (or at the last r-value)
    this->omega = omega_0;
    res = this->integrate(results_0, events, this->get_initial_conditions(), intOpts);
    n_inft_0 = results_0[results_0.size()-1].second[2] > 0._num;    // save if sign(Phi(inf)) is positive or negative

    this->omega = omega_1;
    res_1 = this->integrate(results_1, events, this->get_initial_conditions(), intOpts);
    n_inft_1 = results_1[results_1.size()-1].second[2] > 0._num;
    if (verbose > 0)
        std::cout << "start with omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << " and omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;

    if(res == integrator::endpoint_reached || res_1 == integrator::endpoint_reached) { // in this case psi didn't diverge, so we are probably not integrating to large enough r
        this->r_end *= 1.5_num;
        if (verbose > 1)
            std::cout << "increased r_end to " << this->r_end << " and going deeper" <<  std::endl;
        return bisection(omega_0, omega_1, n_mode, max_steps-steps, delta_omega, verbose);
    }

    steps = 0;
    while((omega_1 - omega_0)/omega_0 > delta_omega && steps < max_steps) { // iterate until accuracy in omega was reached or max number of steps exceeded
        omega_mid = (omega_0 + omega_1)/2._num;
        this->omega = omega_mid;
        res = this->integrate(results_mid, events, this->get_initial_conditions(), intOpts);
        n_inft_mid = results_mid[results_mid.size()-1].second[2] > 0._num;  // save if sign(Phi(inf)) is positive or negative

        if (verbose > 1)
            std::cout << steps << ": omega_mid = " << omega_mid  << " with n_inft= " << n_inft_mid << std::endl;

        // compare the signs of Phi at infinity of the omega-upper, -middle and -lower solution
        // when middle and lower sign are equal, we can move omega_0 to omega_mid
        if(n_inft_mid == n_inft_0) {
            n_inft_0 = n_inft_mid;
            omega_0 = omega_mid;
        }
        else if(n_inft_mid == n_inft_1) {
            n_inft_1 = n_inft_mid;
            omega_1 = omega_mid;
        }
        steps++;
    }

    if (n_inft_1 > 0.)
        this->omega = omega_1;
    else
        this->omega = omega_0;

    if (verbose > 0)
        std::cout << "After " << steps << " steps, chose " << omega << " out of omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << ", omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;

    return 0.;
}

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
int FermionBosonStar::bisection(NUMERIC omega_0, NUMERIC omega_1, int n_mode, int max_steps, NUMERIC delta_omega, int verbose) {

    int result;

    // if phi_0 = 0 then we don't need a bisection
    if (this->phi_0 == 0._num) {
        this->omega = 0._num;
        return 0;
    }
    if(omega_1 < omega_0)
        std::swap(omega_0, omega_1);

    result = this->bisection_find_mode(omega_0, omega_1, n_mode, max_steps, verbose);
    if(result) {
		// failed finding the mode. In some cases this is because the correct mode ends up outside the bisection range during a bisection step
        // re-try the bisection but for different bounds omega_0 and omega_1:
		int max_tries = 20; int tries = 0;
		while (tries < max_tries) {
			if (verbose > -1) {std::cout << "Bisection skipped over the wanted root. Try again with different bounds omega_0 and omega_1. try number:" << tries << std::endl;}
			omega_0 = 1.0; //this->omega * (0.9_num - tries*0.05_num);
			omega_1 = this->omega * (1.1_num + tries*0.31415936_num); // increase the upper range for the bisection by some amount
			result = this->bisection_find_mode(omega_0, omega_1, n_mode, max_steps, verbose);
			tries++;
			if (!result) {if(verbose > -1){std::cout << "re-try was successful after " << tries << " tries" << std::endl;} break;}
		}
	}
	if(result) {return result;}

    result = this->bisection_converge_through_infty_behavior(omega_0, omega_1, n_mode, max_steps, delta_omega, verbose);
    return result;
}


// uses the bisection method to calculate a FBS solution with fixed rho_c and fixed dark matter mass ratio Nb/(Nf+Nb)
// optimized phi_c and omega in the process.
void FermionBosonStar::shooting_NbNf_ratio(NUMERIC NbNf_ratio, NUMERIC NbNf_accuracy, NUMERIC omega_0, NUMERIC omega_1, int n_mode, int max_steps, NUMERIC delta_omega, int verbose) {

	if (NbNf_ratio == 0._num) {
		if (verbose > 0){std::cout << "we have a pure neutron star. No need for shooting" << std::endl;}
		this->phi_0 = 0._num;
		this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
        this->evaluate_model();
		return;
	}
    // before the bisection, calculate the FBS solution once using an initial Phi value
    NUMERIC range_test_for_NbNf;
    NUMERIC phi_c_init = 0.1_num;
	this->phi_0 = phi_c_init;
	if (verbose > 0){std::cout << "start range correction" << std::endl;}
	int i = 0;
    while (i < max_steps) {
        this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
        this->evaluate_model();
        // obtain the ratio Nb/(Nf+Nb)
        range_test_for_NbNf = this->N_B / (this->N_F + this->N_B);
        // check if obtained ratio 'range_test_for_NbNf' is above the wanted ratio:
        // if yes, we perform a bisection search in the range [0, phi_c_init]
        // if no, we increase phi_c_init by an amount until the range covers the wanted value for Nb/(Nf+Nb)
        if (NbNf_ratio < range_test_for_NbNf) {break;}
        else {phi_c_init = phi_c_init*2._num; this->phi_0 = phi_c_init;}
        i++;
        if (verbose > 1) {std::cout << "iteration=" << i << " phi_c_init=" << phi_c_init << " range_test_for_NbNf=" << range_test_for_NbNf << std::endl;}
    }

    // now perform the bisection until the wanted ratio is found with sufficient accuracy:
    NUMERIC phi_c_0 = 0.0_num; 		// lower bound of bisection
    NUMERIC phi_c_1 = phi_c_init; 	// upper bound of bisection
    NUMERIC phi_c_mid = (phi_c_0 + phi_c_1) / 2._num;	// mid point in phi
    NUMERIC NbNf_0 = 0.0_num;	// is zero because scalar field is zero
	NUMERIC NbNf_1 = 100._num;	// initialize with some random large number
    NUMERIC NbNf_mid;

    i = 0;	// reset variable
    // start bisection
	if (verbose > 0){std::cout << "start NbNf-ratio bisection" << std::endl;}
    while ( (std::abs(NbNf_0 - NbNf_1) > NbNf_accuracy) && (i < max_steps) ) {
		
        phi_c_mid = (phi_c_0 + phi_c_1) / 2.;

        this->phi_0 = phi_c_mid;
        this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
        this->evaluate_model();
        NbNf_mid = this->N_B / (this->N_F + this->N_B);

        if (NbNf_mid < NbNf_ratio) {
            		// the mid point is below the wanted ratio and we can adjust the lower bound
            NbNf_0 = NbNf_mid;
            phi_c_0 = phi_c_mid;
			if (verbose > 1){std::cout << "NbNf_mid is smaller than wanted ratio" << std::endl;}
        } else {	// the mid point is above the wanted ratio and we can adjust the upper bound
            NbNf_1 = NbNf_mid;
            phi_c_1 = phi_c_mid;
			if (verbose > 1){std::cout << "NbNf_mid is larger than wanted ratio" << std::endl;}
        }
		i++;
		if (verbose > 1){std::cout << "iteration=" << i << " NbNf_0=" << NbNf_0 << " NbNf_1=" << NbNf_1 << " NbNf_mid=" << NbNf_mid << " NbNf_ratio=" << NbNf_ratio << std::endl;}
    }
}

/* This function takes the result of an integration of the FBS system
 * and calculates the star properties
 * M_T, N_B, N_F, R_B, R_F
 * */
void FermionBosonStar::calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events) {
    const int step_number = results.size();

    /* find the index where the phi field converged (to accurately compute the bosonic radius component later)
     *  */
    bool phi_converged = this->phi_0 <= 0._num;
    int index_phi_converged = 1;

    if (this->phi_0 > 0._num) {
        if(this->R_B_0 >  0._num) { // we artifically set phi to 0 at some point which makes our lifes much easier
            phi_converged = true;
            while(results[index_phi_converged].first < this->R_B_0 && index_phi_converged < step_number-2)
                index_phi_converged++;
        }
        else { // we couldn't successfully set phi to 0 so find the closest thing
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
        }
    }

    // std::cout << "calculate_star_parameters with phi_converged = " << phi_converged << std::endl;

    /*   M_T
     *   Calculate the ADM mass, since we always have convergence, read it from the end
     * M_T := r/2 * ( 1 - 1/(a^2) )
     * */
    NUMERIC M_T = 0._num;
    auto M_func = [&results](int index) { return results[index].first / 2._num * (1._num - 1._num/pow(results[index].second[0], 2)); };
    auto dM_func = [&results, &M_func](int i) { return  (M_func(i+1) - M_func(i))/(results[i+1].first - results[i].first)/ 2._num
                                                        + (M_func(i) - M_func(i-1))/(results[i].first - results[i-1].first)/2._num; };

    if(phi_converged) { // no divergence -> read M_T out at the end
        M_T = M_func(step_number-1);
    }
    else {
        int index_dM_global_minimum = results.size()-3;
        for (int i=results.size()-3; i > index_phi_converged; i--) {
            if(dM_func(index_dM_global_minimum) != dM_func(index_dM_global_minimum)) // NaN prevention
                index_dM_global_minimum = i;
            if(dM_func(i) < dM_func(index_dM_global_minimum)) // find the global minimum
                index_dM_global_minimum = i;
        }

        M_T = M_func(index_dM_global_minimum);
    }

    // std::cout << "min_index_a: " << min_index_a << " min_index_M: " << min_index_dMdr << " min_index_phi: " << min_index_phi << " res_size:" << results.size() << std::endl;

    /*  N_B, N_F
     *  We need to integrate the particle number densities to obtain N_B, N_F */
    std::vector<NUMERIC> r(step_number), N_B_integrand(step_number), N_F_integrand(step_number);
    vector v;
    NUMERIC rho, eps;

    for(unsigned int i = 0; i < results.size(); i++) {
        r[i] = results[i].first;
        v = results[i].second;
        N_B_integrand[i] = 8._num*M_PI * v[0] * this->omega *  v[2] * v[2] * r[i] * r[i] / v[1];  // get bosonic mass (paricle number) for each r
        if (v[4] < P_ns_min || v[4] < this->EOS->min_P())
            rho = 0._num;
        else
            this->EOS->callEOS(rho, eps, v[4]);
        N_F_integrand[i] = 4._num*M_PI * v[0] * rho * r[i] * r[i] ;   // get fermionic mass (paricle number) for each r
    }

    // Integrate
    std::vector<NUMERIC> N_F_integrated, N_B_integrated;
    integrator::cumtrapz(r, N_F_integrand, N_F_integrated);
    integrator::cumtrapz(r, N_B_integrand, N_B_integrated);

    // Find where 99% of N_B,N_F are reached to get the radii
    NUMERIC N_F =  N_F_integrated[step_number-1],
           N_B =  N_B_integrated[index_phi_converged];

    // first find the index in array where 99% is contained
    int i_B = 0, i_F = 0, i_G = 0;
    unsigned int max_index = step_number;
    for(unsigned int i = 1; i < max_index; i++) {
        if(N_B_integrated[i] < 0.99_num*N_B)
            i_B++;
        if(N_F_integrated[i] < 0.99_num*N_F)
            i_F++;
        if(N_F_integrated[i] + N_B_integrated[i] < 0.99_num*(N_B + N_F))
            i_G++;
    }
    // obtain radius from corresponding index
    NUMERIC R_B = r[i_B],
           R_F_0 = r[i_F],
           R_G = r[i_G];

    /*  R_F
     * compute the fermionic radius R_f using the definition where P(R_f)==0
     * iterate the Pressure-array until we find the first point where the pressure is zero */
    NUMERIC R_F = 0.;
    int index_R_F = 0;
    auto P_func = [&results] (int index) {return results[index].second[4];};

    while( P_func(index_R_F) > std::max(P_ns_min, EOS->min_P())  && index_R_F < step_number-1)
        index_R_F++;
    R_F = r[index_R_F];
    N_F = N_F_integrated[index_R_F];

    //std::cout << "M_T = " << M_T << ", N_B = " << N_B << ", R_B = " << R_B << ", N_F = " << N_F << ", R_F = " << R_F << ", R_F_0 = " << R_F_0 << ", N_B/N_F = " << N_B / N_F << std::endl;
    this->M_T = M_T; this->N_B = N_B; this->N_F = N_F; this->R_B = R_B; this->R_F = R_F; this->R_F_0 = R_F_0; this->R_G = R_G;
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

    const bool force_phi_to_zero = false;
    intOpts.save_intermediate = true;
    intOpts.verbose = 0;

    std::vector<integrator::Event> events;

    integrator::Event P_min_reached = FermionBosonStar::P_min_reached; // the presence of this event will increase the accuracy around R_F
    if(this->rho_0 > 0._num)
        events.push_back(P_min_reached);

    this->integrate_and_avoid_phi_divergence(results, events, intOpts, force_phi_to_zero);

    this->calculate_star_parameters(results, events);

    if(!filename.empty()) {
        plotting::save_integration_data(results, {0,1,2,3,4}, {"a", "alpha", "phi", "Psi", "P"}, filename);

        #ifdef DEBUG_PLOTTING
        plotting::plot_evolution(results, events, {0,1,2,3,4}, {"a", "alpha", "phi", "Psi", "P"}, filename.replace(filename.size()-3, 3, "png"));
        matplotlibcpp::legend(); matplotlibcpp::xscale("log"); matplotlibcpp::yscale("log");
        matplotlibcpp::save(filename); matplotlibcpp::close();
        #endif
    }
}

/* Outputs the star parameters and properties */
std::ostream& operator<<(std::ostream& os, const FermionBosonStar& fbs) {
    return os   << fbs.M_T                   << " "   // total gravitational mass
                << fbs.rho_0                 << " "   // central density
                << fbs.phi_0                 << " "   // central scalar field
                << fbs.R_F*1.476625061_num   << " "   // fermionic radius
                << fbs.N_F                   << " "   // number of fermions
                << fbs.R_B*1.476625061_num   << " "   // bosonic radius
                << fbs.R_B_0*1.476625061_num << " "   // phi converged
                << fbs.N_B                   << " "   // number of bosons
                << fbs.N_B / fbs.N_F         << " "   // ratio N_B / N_F
                << fbs.omega                 << " "   // omega
                << fbs.mu                    << " "   // mass mu
                << fbs.lambda                << " "   // self-interaction parameter lambda
                << fbs.R_G*1.476625061_num            // effective gravitational radius
                            ;
}
/* Gives the labels of the values from the output */
std::vector<std::string> FermionBosonStar::labels() {
    return std::vector<std::string> ({"M_T", "rho_0", "phi_0", "R_F", "N_F", "R_B", "R_B_0", "N_B", "N_B/N_F", "omega", "mu", "lambda", "R_G"});
}
