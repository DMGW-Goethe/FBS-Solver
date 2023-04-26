#include "fps.hpp"


/*   Events   */
/* Overwrite events from the parent class so that we can use the same implementation: */
/* This event triggers when E or B are diverging, i.e. E,B > 1. */
const integrator::Event FermionProcaStar::Psi_diverging = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return (std::abs(y[2]) > 1.0 || std::abs(y[3]) > 1.0); }, true, "E_B_diverging");

/* This event triggers when E becomes negative, it is used to count the zero crossings of E */
const integrator::Event FermionProcaStar::phi_negative = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return y[2] < 0.; }, false, "E_negative");
/* This event triggers when E becomes positive */
const integrator::Event FermionProcaStar::phi_positive = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                                                { return y[2] > 0.; }, false, "E_positive");

/* This event checks that E and B have converged */
const integrator::Event FermionProcaStar::EB_converged = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params)
                                                                            { return (	std::abs(y[2])/((FermionProcaStar*)params)->E_0 < PHI_converged &&
                                                                                    	std::abs(y[3])/((FermionProcaStar*)params)->E_0/*/r*/ < PHI_converged ); }, false, "EB_converged");



/* This function gives the system of ODEs for the FPS star
 *  for the variables a, alpha, E, B, and P
 *  taken from [...]
 *
 *  This function is called by the integrator during the integration
 * */
vector FermionProcaStar::dy_dr(const double r, const vector& vars) const {

    // rename input & class variables for simpler use
    const double a = vars[0]; const double alpha = vars[1]; const double E = vars[2]; const double B = vars[3]; double P = vars[4];
    EquationOfState& EoS = *(this->EOS);

    // define hydrodynamic quantities
    double etot = 0.;	// total energy density. Must be computed from P using the EoS
    // define potential of the bosonic field
	const double AmuAmu = B*B/a/a - E*E/alpha/alpha;	// contraction term: A_mu*A^mu
    const double V = mu*mu*AmuAmu + lambda/2.*AmuAmu*AmuAmu;
    const double dV_deps = mu*mu + lambda*AmuAmu;
	const double dV2_deps2 = lambda;

    // apply the EOS
    if(P <= 0. || P < EoS.min_P() || P < P_ns_min)  {
        P = 0.; etot = 0.;
    } else {
        etot = EoS.get_e_from_P(P);
    }

    // compute the ODEs:
	#define DEF_KAPPA 1.0
	//8.*M_PI
	//DEF_KAPPA is there to use units where G=1/8pi
	double dE_dr = - dV_deps*B*alpha*alpha/omega + omega*B;
    double da_dr = 0.5* a *      ( (1.-a*a) / r + DEF_KAPPA*r*a*a*( etot + std::pow((dE_dr - omega*B),2)/(alpha*alpha*a*a) + V + 2*dV_deps*E*E/alpha/alpha ));
    double dalpha_dr = 0.5* alpha * ( (a*a-1.) / r + DEF_KAPPA*r*a*a*( P - std::pow((dE_dr - omega*B),2)/(alpha*alpha*a*a) - V + 2*dV_deps*B*B/a/a) );
    
	double dB_dr = ( dV2_deps2*( 2*B*B*da_dr/a/a/a + 2.*E*dE_dr/alpha/alpha - 2.*E*E*dalpha_dr/alpha/alpha/alpha)*B*alpha*alpha  
					- dV_deps*(a*a*E*omega + 2.*B*alpha*dalpha_dr)
					- ( da_dr/a + dalpha_dr/alpha - 2./r)*omega*(dE_dr - omega*B) )
					/ (dV2_deps2*2.*B*B*alpha*alpha/(a*a) + dV_deps*alpha*alpha);
    
	double dP_dr = -(etot + P)*dalpha_dr/alpha;

    // write the ODE values into output vector
    return vector({da_dr, dalpha_dr, dE_dr, dB_dr, dP_dr});
}


/* This function gives the initial conditions for the FBS integration
 * with a_0 = 1,  alpha_0 = 1,  E_0 = this->E_0, B = 0, and P_0 = EOS(rho_0)
 * at the position r_init (which is assumed to be small enough for this to be valid) */
vector FermionProcaStar::get_initial_conditions(double r_init) const {
    return vector( {1.0, 1.0, this->E_0, 0., rho_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(this->rho_0, 0.) : 0.});
}

/* FUNCTIONS RELATED TO THE BISECTION, TO FIND THE RIGHT MODE */

/* This function performs an integration and tries to find R_B_0 while integrating, where the bosonic component has sufficiently converged
 *  and we can set it to zero. The convergence criterion is given by the EB_converged event. The function returns R_B_0 by reference.
 * If the condition is not fulfilled, the "convergence" can be forced (with the force boolean) by finding the last minimum of the phi field.
 * */
int FermionProcaStar::find_bosonic_convergence(std::vector<integrator::step>& results, std::vector<integrator::Event>& events, integrator::IntegrationOptions intOpts, double& R_B_0, bool force, double r_init, double r_end) const {

    if(this->phi_0 <= 0.)
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
            y_at_EB_converged[3] =0.;
            if (!EB_converged.condition(results[index_phi_converged].first, 0., y_at_EB_converged, dy_at_phi_converged, (const void*)this)) {
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


/* FUNCTIONS TO COMPUTE MACROSCOPIC PARAMETERS OF THE STAR */

/* This function takes the result of an integration of the FPS system
 * and calculates the star properties
 * M_T, N_B, N_F, R_B, R_F, etc.
 * */
void FermionProcaStar::calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events) {
    const int step_number = results.size();

    /* find the index where the phi field converged (to accurately compute the bosonic radius component later)
     *  */
    bool E_converged = this->E_0 <= 0.;
    int index_E_converged = 1;

    if (this->phi_0 > 0.) {
        if(this->R_B_0 >  0.) { // we artifically set phi to 0 at some point which makes our lifes much easier
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
    double M_T = 0.;
    auto M_func = [&results](int index) { return results[index].first / 2. * (1. - 1./pow(results[index].second[0], 2)); };
    auto dM_func = [&results, &M_func](int i) { return  (M_func(i+1) - M_func(i))/(results[i+1].first - results[i].first)/ 2.
                                                        + (M_func(i) - M_func(i-1))/(results[i].first - results[i-1].first)/2.; };

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
    std::vector<double> r(step_number), N_B_integrand(step_number), N_F_integrand(step_number);
    vector v;
    double rho, eps;

    for(unsigned int i = 0; i < results.size(); i++) {
        r[i] = results[i].first;
        v = results[i].second;
		double dV = this->mu*this->mu + this->lambda * (v[3]*v[3]/v[0]/v[0] - v[2]*v[2]/v[1]/v[1]);
		double dE_dr = this->omega* v[3] - dV * v[3] * v[1]*v[1] / this->omega;
        N_B_integrand[i] = 8.*M_PI * v[3] * (this->omega * v[3] - dE_dr) * r[i] * r[i] / (v[0]*v[1]);  // get bosonic mass (paricle number) for each r
        if (v[4] < P_ns_min || v[4] < this->EOS->min_P())
            rho = 0.;
        else
            this->EOS->callEOS(rho, eps, v[4]);
        N_F_integrand[i] = 4.*M_PI * v[0] * rho * r[i] * r[i];   // get fermionic mass (paricle number) for each r
    }

    // Integrate
    std::vector<double> N_F_integrated, N_B_integrated;
    integrator::cumtrapz(r, N_F_integrand, N_F_integrated);
    integrator::cumtrapz(r, N_B_integrand, N_B_integrated);

    // Find where 99% of N_B,N_F are reached to get the radii
    double N_F =  N_F_integrated[step_number-1],
           N_B =  N_B_integrated[index_E_converged];

    // first find the index in array where 99% is contained
    int i_B = 0, i_F = 0, i_G = 0;
    unsigned int max_index = step_number;
    for(unsigned int i = 1; i < max_index; i++) {
        if(N_B_integrated[i] < 0.99*N_B)
            i_B++;
        if(N_F_integrated[i] < 0.99*N_F)
            i_F++;
        if(N_F_integrated[i] + N_B_integrated[i] < 0.99*(N_B + N_F))
            i_G++;
    }
    // obtain radius from corresponding index
    double R_B = r[i_B],
           R_F_0 = r[i_F],
           R_G = r[i_G];

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
    if(this->rho_0 > 0.)
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
                << fps.R_F*1.476625061       << " "   // fermionic radius
                << fps.N_F                   << " "   // number of fermions
                << fps.R_B*1.476625061       << " "   // bosonic radius
                << fps.R_B_0*1.476625061     << " "   // phi converged
                << fps.N_B                   << " "   // number of bosons
                << fps.N_B / fps.N_F         << " "   // ratio N_B / N_F
                << fps.omega                 << " "   // omega
                << fps.mu                    << " "   // mass mu
                << fps.lambda                << " "   // self-interaction parameter lambda
                << fps.R_G*1.476625061                // effective gravitational radius
                            ;
}
/* Gives the labels of the values from the output */
std::vector<std::string> FermionProcaStar::labels() {
    return std::vector<std::string> ({"M_T", "rho_0", "E_0", "R_F_0", "N_F", "R_B", "R_B_0", "N_B", "N_B/N_F", "omega", "mu", "lambda", "R_G"});
}
