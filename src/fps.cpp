#include "fps.hpp"


/*   Events   */
/* Overwrite events from the parent class so that we can use the same implementation: */
/* This event triggers when E or B are diverging, i.e. E,B > 1. */
const integrator::Event FermionProcaStar::Psi_diverging = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (std::abs(y[2]) > 1.0_num || std::abs(y[3]) > 1.0_num); }, true, "E_B_diverging");

/* This event triggers when E becomes negative, it is used to count the zero crossings of E */
const integrator::Event FermionProcaStar::phi_negative = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return y[2] < 0._num; }, false, "E_negative");
/* This event triggers when E becomes positive */
const integrator::Event FermionProcaStar::phi_positive = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return y[2] > 0._num; }, false, "E_positive");

/* This event checks that E and B have converged */
const integrator::Event FermionProcaStar::EB_converged = integrator::Event([](const NUMERIC r, const NUMERIC dr, const vector& y, const vector& dy, const void*params)
                                                                            { return (	std::abs(y[2])/((FermionProcaStar*)params)->E_0 < PHI_converged &&
                                                                                    	std::abs(y[3])/((FermionProcaStar*)params)->E_0/*/r*/ < PHI_converged ); }, false, "EB_converged");



/* This function gives the system of ODEs for the FPS star
 *  for the variables a, alpha, E, B, and P
 *  taken from [...]
 *
 *  This function is called by the integrator during the integration
 * */
vector FermionProcaStar::dy_dr(const NUMERIC r, const vector& vars) const {

    // rename input & class variables for simpler use
    const NUMERIC a = vars[0]; const NUMERIC alpha = vars[1]; const NUMERIC E = vars[2]; const NUMERIC B = vars[3]; NUMERIC P = vars[4];
    EquationOfState& EoS = *(this->EOS);

    // define hydrodynamic quantities
    NUMERIC etot = 0._num;	// total energy density. Must be computed from P using the EoS
    // define potential of the bosonic field
	const NUMERIC AmuAmu = B*B/a/a - E*E/alpha/alpha;	// contraction term: A_mu*A^mu
    const NUMERIC V = mu*mu*AmuAmu + lambda/2.*AmuAmu*AmuAmu;
    const NUMERIC dV_deps = mu*mu + lambda*AmuAmu;
	const NUMERIC dV2_deps2 = lambda;

    // apply the EOS
    if(P <= 0._num || P < EoS.min_P() || P < P_ns_min)  {
        P = 0._num; etot = 0._num;
    } else {
        etot = EoS.get_e_from_P(P);
    }

    // compute the ODEs:
	#define DEF_KAPPA 8._num*M_PI
	//8.*M_PI
	//DEF_KAPPA is there to use units where G=1/8pi
	NUMERIC dE_dr = - dV_deps*B*alpha*alpha/omega + omega*B;
	if (omega == 0._num) {dE_dr=0._num;}	// need this failsafe when computing pure NS
    NUMERIC da_dr = 0.5_num* a *      ( (1._num-a*a) / r + DEF_KAPPA*r*a*a*( etot + std::pow((dE_dr - omega*B),2)/(alpha*alpha*a*a) + V + 2._num*dV_deps*E*E/alpha/alpha ));
    NUMERIC dalpha_dr = 0.5_num* alpha * ( (a*a-1._num) / r + DEF_KAPPA*r*a*a*( P - std::pow((dE_dr - omega*B),2)/(alpha*alpha*a*a) - V + 2._num*dV_deps*B*B/a/a) );
    
	NUMERIC dB_dr = ( dV2_deps2*( 2._num*B*B*da_dr/a/a/a + 2._num*E*dE_dr/alpha/alpha - 2._num*E*E*dalpha_dr/alpha/alpha/alpha)*B*alpha*alpha  
					- dV_deps*(a*a*E*omega + 2._num*B*alpha*dalpha_dr)
					- ( da_dr/a + dalpha_dr/alpha - 2._num/r)*omega*(dE_dr - omega*B) )
					/ (dV2_deps2*2._num*B*B*alpha*alpha/(a*a) + dV_deps*alpha*alpha);
    
	NUMERIC dP_dr = -(etot + P)*dalpha_dr/alpha;

    // write the ODE values into output vector
    return vector({da_dr, dalpha_dr, dE_dr, dB_dr, dP_dr});
}


/* This function gives the initial conditions for the FBS integration
 * with a_0 = 1,  alpha_0 = 1,  E_0 = this->E_0, B = 0, and P_0 = EOS(rho_0)
 * at the position r_init (which is assumed to be small enough for this to be valid) */
vector FermionProcaStar::get_initial_conditions(NUMERIC r_init) const {
    return vector( {1._num, 1._num, this->E_0, 0._num, rho_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(this->rho_0, 0._num) : 0._num});
}

/* FUNCTIONS RELATED TO THE BISECTION, TO FIND THE RIGHT MODE */

/* This function performs an integration and tries to find R_B_0 while integrating, where the bosonic component has sufficiently converged
 *  and we can set it to zero. The convergence criterion is given by the EB_converged event. The function returns R_B_0 by reference.
 * If the condition is not fulfilled, the "convergence" can be forced (with the force boolean) by finding the last minimum of the phi field.
 * */
int FermionProcaStar::find_bosonic_convergence(std::vector<integrator::step>& results, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts, NUMERIC& R_B_0, bool force, NUMERIC r_init, NUMERIC r_end) const {

    if(this->phi_0 <= 0._num)
        return -1;

    //integrator::Event Psi_positive = FermionBosonStar::Psi_positive;
    integrator::Event EB_converged = FermionProcaStar::EB_converged;
    EB_converged.stopping_condition = true;
    //events.push_back(Psi_positive); // the presence of this event increases the accuracy around zero crossings
    events.push_back(EB_converged);

    int res = this->integrate(results, events, this->get_initial_conditions(), intOpts, r_init, r_end);

    EB_converged = events[events.size()-1];
    events.pop_back(); // Remove events that we added
    //events.pop_back();

    if (res == integrator::event_stopping_condition && !EB_converged.active) // in this case another event triggered the stop so just return
        return res;

    if (!EB_converged.active ) { // the integrator didn't catch the convergence so we have to find it ourselves

        if (!intOpts.save_intermediate) { // can't find it if we don't have the steps
            results.clear();
            intOpts.save_intermediate = true;
            res = this->integrate(results, events, this->get_initial_conditions(), intOpts, r_init, r_end);
            intOpts.save_intermediate = false;
        }
        // find the minimum in the phi field before it diverges
        int index_phi_converged = 1;

        auto abs_EB_func = [&results] (int index) { return std::abs(results[index].second[2]) + std::abs(results[index].second[3]); };
		std::vector<int> abs_phi_minima({});
        int index_phi_global_min = 0;
        for (unsigned int i=1; i < results.size()-1; i++) {
            if ( abs_EB_func(i) <= abs_EB_func(i-1) && abs_EB_func(i) < abs_EB_func(i+1) )
                abs_phi_minima.push_back(i);
            if (abs_EB_func(i) < abs_EB_func(index_phi_global_min) )
                index_phi_global_min = i;
        }
        index_phi_converged = index_phi_global_min;
        if(abs_phi_minima.size() > 0)
            index_phi_converged = abs_phi_minima[abs_phi_minima.size()-1];

        // maybe the event didn't trigger because psi was too large?
        if(!force) {
            auto y_at_EB_converged = results[index_phi_converged].second;
            vector dy_at_phi_converged = this->dy_dr(results[index_phi_converged].first, y_at_EB_converged);
            y_at_EB_converged[3] =0._num;
            if (!EB_converged.condition(results[index_phi_converged].first, 0._num, y_at_EB_converged, dy_at_phi_converged, (const void*)this)) {
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

int FermionProcaStar::integrate_and_avoid_phi_divergence(std::vector<integrator::step>& results, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts, bool force, std::vector<int> additional_zero_indices, NUMERIC r_init, NUMERIC r_end)  {

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


// uses the bisection method to calculate a FBS solution with fixed rho_c and fixed dark matter mass ratio Nb/(Nf+Nb)
// optimized phi_c and omega in the process.
void FermionProcaStar::shooting_NbNf_ratio(NUMERIC NbNf_ratio, NUMERIC NbNf_accuracy, NUMERIC omega_0, NUMERIC omega_1, int n_mode, int max_steps, NUMERIC delta_omega, int verbose) {

	if (NbNf_ratio == 0._num) {
		if (verbose > 0){std::cout << "we have a pure neutron star. No need for shooting" << std::endl;}
		this->E_0 = 0._num; this->phi_0 = this->E_0;
		this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
        this->evaluate_model();
		return;
	}
    // before the bisection, find the ideal upper and lower bounds for E_0. Also chek for edge cases:
	// beware: E_0 MUST have: E_0 < m / sqrt(lambda) = 1/4*sqrt(pi*Lambda_int)
    NUMERIC range_test_for_NbNf;
    NUMERIC E0_init = 0.1_num;

	if (verbose > 0){std::cout << "start upper bound correction for E_0" << std::endl;}
	int i = 0;
	if (abs(lambda) > 0._num) {	// a nonzero self-interaction produces an upper bound for the central value of E_0
		E0_init = 0.999_num*mu/std::sqrt(lambda);	// maximum allowed value for E_0
		this->E_0 = E0_init; this->phi_0 = this->E_0;
		// chek if the maximum possible value of E_0 is able to produce a high enough Nb/(Nf+Nb)-ratio for the bisection:
		// if yes, we perform a bisection search in the range [0, E0_init]
		// if no, we cannot do the bisection and exit this function
		this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
		this->evaluate_model();
		range_test_for_NbNf = this->N_B / (this->N_F + this->N_B);
		if (NbNf_ratio > range_test_for_NbNf) {
			if (verbose > 0){std::cout << "maximal possible E_0 is not large enough to produce the wanted Nb/(Nf+Nb)-ratio" << std::endl;}
			this->E_0 = 0._num; this->phi_0 = 0._num; this->rho_0 = 0._num; // we obtained a solution with not the wanted ratio. Set star params to zero to mark the wrong results
			return;
		}
	} 
	else {	// self-interaction is zero, there is no upper bound for E_0 and we can chose a suitable value
		this->E_0 = E0_init; this->phi_0 = this->E_0;
		if (verbose > 0){std::cout << "start range correction" << std::endl;}
    	while (i < max_steps) {
			this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
			this->evaluate_model();
			// obtain the ratio Nb/(Nf+Nb)
			range_test_for_NbNf = this->N_B / (this->N_F + this->N_B);
			// check if obtained ratio 'range_test_for_NbNf' is above the wanted ratio:
			// if yes, we perform a bisection search in the range [0, E0_init]
			// if no, we increase E0_init by an amount until the range covers the wanted value for Nb/(Nf+Nb)
			if (NbNf_ratio < range_test_for_NbNf) {break;}
			else {E0_init = E0_init*2._num; this->E_0 = E0_init; this->phi_0 = this->E_0;}
			i++;
			if (verbose > 1) {std::cout << "iteration=" << i << " E0_init=" << E0_init << " range_test_for_NbNf=" << range_test_for_NbNf << std::endl;}
    	}
	}
    // now perform the bisection until the wanted ratio is found with sufficient accuracy:
    NUMERIC E0_0 = 0.0_num; 		// lower bound of bisection
    NUMERIC E0_1 = E0_init; 	// upper bound of bisection
    NUMERIC E0_mid = (E0_0 + E0_1) / 2._num;	// mid point in phi
    NUMERIC NbNf_0 = 0.0_num;	// is zero because scalar field is zero
	NUMERIC NbNf_1 = 100._num;	// initialize with some random large number
    NUMERIC NbNf_mid;

    i = 0;	// reset variable
    // start bisection
	if (verbose > 0){std::cout << "start NbNf-ratio bisection" << std::endl;}
    while ( (std::abs(NbNf_0 - NbNf_1) > NbNf_accuracy) && (i < max_steps) ) {
		
        E0_mid = (E0_0 + E0_1) / 2.;

        this->E_0 = E0_mid; this->phi_0 = this->E_0;
        this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
        this->evaluate_model();
        NbNf_mid = this->N_B / (this->N_F + this->N_B);

        if (NbNf_mid < NbNf_ratio) {
            		// the mid point is below the wanted ratio and we can adjust the lower bound
            NbNf_0 = NbNf_mid;
            E0_0 = E0_mid;
			if (verbose > 1){std::cout << "NbNf_mid is smaller than wanted ratio" << std::endl;}
        } else {	// the mid point is above the wanted ratio and we can adjust the upper bound
            NbNf_1 = NbNf_mid;
            E0_1 = E0_mid;
			if (verbose > 1){std::cout << "NbNf_mid is larger than wanted ratio" << std::endl;}
        }
		i++;
		if (verbose > 1){std::cout << "iteration=" << i << " NbNf_0=" << NbNf_0 << " NbNf_1=" << NbNf_1 << " NbNf_mid=" << NbNf_mid << " NbNf_ratio=" << NbNf_ratio << std::endl;}
    }
}

/* FUNCTIONS TO COMPUTE MACROSCOPIC PARAMETERS OF THE STAR */

/* This function takes the result of an integration of the FPS system
 * and calculates the star properties
 * M_T, N_B, N_F, R_B, R_F, etc.
 * */
void FermionProcaStar::calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events) {
    const int step_number = results.size();

    /* find the index where the phi field converged (to accurately compute the bosonic radius component later)
     *  */
    bool E_converged = this->E_0 <= 0._num;
    int index_E_converged = 1;

    if (this->phi_0 > 0._num) {
        if(this->R_B_0 >  0._num) { // we artifically set phi to 0 at some point which makes our lifes much easier
            E_converged = true;
            while(results[index_E_converged].first < this->R_B_0 && index_E_converged < step_number-2)
                index_E_converged++;
        }
        else { // we couldn't successfully set E to 0 so find the closest thing
            auto abs_E_func = [&results] (int index) { return std::abs(results[index].second[2]); };
            std::vector<int> abs_E_minima({});
            int index_E_global_min = 0;
            for (unsigned int i=1; i < results.size()-1; i++) {
                if ( abs_E_func(i) <= abs_E_func(i-1) && abs_E_func(i) < abs_E_func(i+1) )
                    abs_E_minima.push_back(i);
                if (abs_E_func(i) < abs_E_func(index_E_global_min) )
                    index_E_global_min = i;
            }
            index_E_converged = index_E_global_min;
            if(abs_E_minima.size() > 0)
                index_E_converged = abs_E_minima[abs_E_minima.size()-1];
        }
    }

    // std::cout << "calculate_star_parameters with E_converged = " << E_converged << std::endl;

    /*   M_T
     *   Calculate the ADM mass, since we always have convergence, read it from the end
     * M_T := r/2 * ( 1 - 1/(a^2) )
     * */
    NUMERIC M_T = 0._num;
    auto M_func = [&results](int index) { return results[index].first / 2._num * (1._num - 1._num/pow(results[index].second[0], 2)); };
    auto dM_func = [&results, &M_func](int i) { return  (M_func(i+1) - M_func(i))/(results[i+1].first - results[i].first)/ 2._num
                                                        + (M_func(i) - M_func(i-1))/(results[i].first - results[i-1].first)/2._num; };

    if(E_converged) { // no divergence -> read M_T out at the end
        M_T = M_func(step_number-1);
    }
    else {
        int index_dM_global_minimum = results.size()-3;
        for (int i=results.size()-3; i > index_E_converged; i--) {
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
		NUMERIC dV = this->mu*this->mu + this->lambda * (v[3]*v[3]/v[0]/v[0] - v[2]*v[2]/v[1]/v[1]);
		NUMERIC dE_dr = this->omega* v[3] - dV * v[3] * v[1]*v[1] / this->omega;
        N_B_integrand[i] = 8._num*M_PI * v[3] * (this->omega * v[3] - dE_dr) * r[i] * r[i] / (v[0]*v[1]);  // get bosonic mass (paricle number) for each r
        if (v[4] < P_ns_min || v[4] < this->EOS->min_P())
            rho = 0._num;
        else
            this->EOS->callEOS(rho, eps, v[4]);
        N_F_integrand[i] = 4._num*M_PI * v[0] * rho * r[i] * r[i];   // get fermionic mass (paricle number) for each r
    }

    // Integrate
    std::vector<NUMERIC> N_F_integrated, N_B_integrated;
    integrator::cumtrapz(r, N_F_integrand, N_F_integrated);
    integrator::cumtrapz(r, N_B_integrand, N_B_integrated);

    // Find where 99% of N_B,N_F are reached to get the radii
    NUMERIC N_F =  N_F_integrated[step_number-1],
           N_B =  N_B_integrated[index_E_converged];

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
    NUMERIC R_F = 0._num;
    int index_R_F = 0;
    auto P_func = [&results] (int index) {return results[index].second[4];};

    while( P_func(index_R_F) > std::max(P_ns_min, EOS->min_P())  && index_R_F < step_number-1)
        index_R_F++;
    R_F = r[index_R_F];
    N_F = N_F_integrated[index_R_F];

    //std::cout << "M_T = " << M_T << ", N_B = " << N_B << ", R_B = " << R_B << ", N_F = " << N_F << ", R_F = " << R_F << ", R_F_0 = " << R_F_0 << ", N_B/N_F = " << N_B / N_F << std::endl;
    this->M_T = M_T; this->N_B = N_B; this->N_F = N_F; this->R_B = R_B; this->R_F = R_F; this->R_F_0 = R_F_0; this->R_G = R_G;
}


/* TOP-LEVEL LOGIC: FUNCTIONS RELATED TO THE FINAL EVALUATION OF THE SOLUTION*/

/* Simple wrapper function if the results of the integration are not needed for the user */
void FermionProcaStar::evaluate_model() {
    std::vector<integrator::step> results;
    this->evaluate_model(results);
}

/* This function integrates over the ODE system while avoiding the phi divergence
 *  and then calculates the star properties
 *  Optionally, the output of the integration is saved in the file
 *  Only call if omega is the corresponding eigenfrequency of mu, lambda, rho_0, phi_0 */
void FermionProcaStar::evaluate_model(std::vector<integrator::step>& results, integrator::IntegrationOptions intOpts, std::string filename) {

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
        plotting::save_integration_data(results, {0,1,2,3,4}, {"a", "alpha", "E", "B", "P"}, filename);
    }
}

/* FUNCTIONS RELATED TO OUTPUT OF VARIABLES */

/* Outputs the star parameters and properties */
std::ostream& operator<<(std::ostream& os, const FermionProcaStar& fps) {
    return os   << fps.M_T                   << " "   // total gravitational mass
                << fps.rho_0                 << " "   // central density
                << fps.E_0                   << " "   // central vector field A_t component
                << fps.R_F*1.476625061_num   << " "   // fermionic radius
                << fps.N_F                   << " "   // number of fermions
                << fps.R_B*1.476625061_num   << " "   // bosonic radius
                << fps.R_B_0*1.476625061_num << " "   // phi converged
                << fps.N_B                   << " "   // number of bosons
                << fps.N_B / fps.N_F         << " "   // ratio N_B / N_F
                << fps.omega                 << " "   // omega
                << fps.mu                    << " "   // mass mu
                << fps.lambda                << " "   // self-interaction parameter lambda
                << fps.R_G*1.476625061_num            // effective gravitational radius
                            ;
}
/* Gives the labels of the values from the output */
std::vector<std::string> FermionProcaStar::labels() {
    return std::vector<std::string> ({"M_T", "rho_0", "E_0", "R_F_0", "N_F", "R_B", "R_B_0", "N_B", "N_B/N_F", "omega", "mu", "lambda", "R_G"});
}
