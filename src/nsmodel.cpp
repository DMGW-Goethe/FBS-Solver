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

    if(P < 0.) P = 1e-20;  // need this to prevent NaN errors...

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


void FermionBosonStar::set_initial_conditions(const double rho_0, const double phi_0) {
    this->rho_0 = rho_0;
    this->phi_0 = phi_0;
    this->initial_conditions =  vector( {1.0, 1.0, phi_0, 0., this->EOS->get_P_from_rho(rho_0)});
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
    double r_init = 1e-10, r_end= 1000;     // initial integration radius and max radius for integration
    // define events to check for during integration and put them inside of a std::vector:
    integrator::Event phi_neg([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[2] < 0.; });
    integrator::Event phi_pos([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[2] > 0.; });
    // stop integration if solution diverges:
    integrator::Event sol_diverged([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return std::abs(y[3]) > 1.0; }, true);
    std::vector<integrator::Event> events = {phi_neg, phi_pos, sol_diverged};     // put the events into the event array
    // declare containers to hold the solution of the integration for the upper- (1), lower- (0) and middle (mid) omega
    std::vector<integrator::step> results_0, results_1, results_mid;

    // find initial values for omega min and omega max
    assert(omega_0 < omega_1);  // if the lower omega (omega_0) is larger than the upper omega (omega_1), we stop the code execution

    // set the lower omega and integrate the ODEs:
    this->omega = omega_0;
    int res =  integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  results_0,  events, intOpts);
    n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    // set the upper omega and integrate the ODEs:
    this->omega = omega_1;
    res =  integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  results_1,  events, intOpts);
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
        res =  integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  results_mid,  events, intOpts);
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
    res =  integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  results_0,  events, intOpts);
    n_inft_0 = results_0[results_0.size()-1].second[2] > 0.;    // save if sign(Phi(inf)) is positive or negative

    this->omega = omega_1;
    res =  integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  results_1,  events, intOpts);
    n_inft_1 = results_1[results_1.size()-1].second[2] > 0.;    // save if sign(Phi(inf)) is positive or negative
    //std::cout << "start with omega_0 =" << omega_0 << " with n_inft=" << n_inft_0 << " and omega_1=" << omega_1 << " with n_inft=" << n_inft_1 << std::endl;

    /* plot the evolution with python
    #ifdef DEBUG_PLOTTING
    plotting::plot_evolution(results_0, events, {2}, { "Phi_0" });
    plotting::plot_evolution(results_1, events, {2}, { "Phi_1" });
    matplotlibcpp::legend();
    matplotlibcpp::save("bisection_int.png"); matplotlibcpp::close();
    #endif
    */

    intOpts.save_intermediate=false;
    while(omega_1 - omega_0 > delta_omega && i < max_steps) { // iterate until accuracy in omega was reached or max number of steps exceeded
        omega_mid = (omega_0 + omega_1)/2.;
        //std::cout << "omega_mid = " << omega_mid << " ->";
        this->omega = omega_mid;
        res =  integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  results_mid,  events, intOpts);
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

    /* integrate again with saving the intermediate steps to plot them (we remove this later)
    this->omega = omega_0; intOpts.save_intermediate=true;
    res =  integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  results_0,  events, intOpts);

    this->omega = omega_1;
    res =  integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  results_1,  events, intOpts);

    #ifdef DEBUG_PLOTTING
    plotting::plot_evolution(results_0, events, {2}, {"Phi_0"});
    plotting::plot_evolution(results_1, events, {2}, {"Phi_1"});
    matplotlibcpp::legend();
    matplotlibcpp::save("bisection_fin.png"); matplotlibcpp::close();
    #endif
    */

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

        this->initial_conditions[2] = phi_c_init;
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

    this->initial_conditions[2] = phi_c_0;
    this->bisection(omega_0, omega_1, n_mode, max_steps, delta_omega);
    this->evaluate_model();
    NbNf_0 = this->N_B / this->N_F;

    int i = 0;
    // continue bisection until the wanted accuracy was reached
    while ( (std::abs(NbNf_0 - NbNf_1) > NbNf_accuracy) && (i < max_steps) ) {
        i++;

        phi_c_mid = (phi_c_0 + phi_c_1) / 2.;

        this->initial_conditions[2] = phi_c_mid;
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
    double r_init = 1e-10, r_end= 1000;

    integrator::Event M_converged([](const double r, const double dr, const vector& y, const vector& dy, const void *params) {
                                                                                                        double dM_dr = ((1. - 1./y[0]/y[0])/2. + r*dy[0]/y[0]/y[0]/y[0]);
                                                                                                        return  dM_dr < 1e-18 ; },true);
    // stop integration if solution diverges:
    integrator::Event sol_diverged([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return (std::abs(y[3]) > 1.0); }, true);

    std::vector<integrator::Event> events = {M_converged, sol_diverged};
    //std::vector<integrator::step> results;
    results.clear();

    int res =  integrator::RKF45(&(this->dy_dt_static), r_init, this->initial_conditions, r_end, (void*) this,  results,  events, intOpts);

    if(!filename.empty()) {
        plotting::save_integration_data(results, {0,1,2,3,4}, {"a", "alpha", "Phi", "Psi", "P"}, filename);

        #ifdef DEBUG_PLOTTING
        plotting::plot_evolution(results, events, {0,1,2,3,4}, {"a", "alpha", "Phi", "Psi", "P"}, filename.replace(filename.size()-3, 3, "png"));
        matplotlibcpp::legend(); matplotlibcpp::yscale("log");
        matplotlibcpp::save(filename); matplotlibcpp::close();
        #endif
    }

    // obtain estimate for the total mass:
    // find the minimum in the g_tt component and then compute the total mass M_T at the corresponding index:
    // start iterating through the solution array backwards
    int min_index_a = results.size()-1;
    double curr_a_min = results[results.size()-1].second[0];  // last value for the metric component a at r=0
    for (int i=results.size()-2; i >= 0; i--) {
        if (results[i].second[0] < curr_a_min) {
            curr_a_min = results[i].second[0]; // update current minimum
            min_index_a = i;  // update min_index
        }
        else {
            break; // the component is increasing again. that means we have found out minimum
        }
    }
    // find the minimum in the phi field before it diverges (to accurately compute the bosonic radius component later):
    // Note: this method will maybe not work well if we consider higher modes of phi!
    int min_index_phi = results.size()-1;
    double curr_phi_min = std::abs(results[results.size()-1].second[2]);  // last value for the phi field a at r=0
    for (int i=results.size()-2; i >= 0; i--) {
        if (std::abs(results[i].second[2]) < curr_phi_min) {
            curr_phi_min = std::abs(results[i].second[2]); // update current minimum
            min_index_phi = i;  // update min_index
        }
        else {
            break; // the component is increasing again. that means we have found out minimum
        }
    }

    // calculate M_T in where the minimum of the 'a' metric component is:
    // M_T = r/2 * (1 - 1/a^2)
    double M_T = results[min_index_a].first / 2. * (1. - 1./results[min_index_a].second[0]/results[min_index_a].second[0]);

    std::cout << "min_index_a: " << min_index_a << " min_index_phi: " << min_index_phi << " res_size:" << results.size() << std::endl;

    // Extract the results and put them into a usable form to calculate N_B, N_F
    std::vector<double> r(results.size()), N_B_integrand(results.size()), N_F_integrand(results.size());
    vector v;
    double rho, eps;

    for(int i = 0; i < results.size(); i++) {
        r[i] = results[i].first;
        v = results[i].second;
        N_B_integrand[i] = v[0] * this->omega *  v[2] * v[2] * r[i] * r[i] / v[1];  // get bosonic mass (paricle number) for each r
        this->EOS->callEOS(rho, eps, std::max(0., v[4]));
        N_F_integrand[i] = v[0] * rho * r[i] * r[i] ;   // get fermionic mass (paricle number) for each r
    }

    // Integrate
    std::vector<double> N_F_integrated, N_B_integrated;
    integrator::cumtrapz(r, N_F_integrand, N_F_integrated);
    integrator::cumtrapz(r, N_B_integrand, N_B_integrated);

    // Find where 99% of N_B,N_F are reached to get the radii
    // we must take the value of the integral *before* the solution diverges! Therefore we cannot just take the last array element
    // but we take the index of the minimum of the metrig g_tt component and the scalar field Phi respectively! This is given by "min_index_*" (see above)
    double N_F = N_F_integrated[min_index_a], N_B = N_B_integrated[min_index_phi];

    // first find the index in array where 99% is contained
    // only iterate until the position where the minimum of the metrig g_tt component is (min_index)
    int i_B=-1, i_F=-1; // index of fermionic/bsonic radius
    for(int i = 0; i < results.size(); i++) {
        if(i_B < 0) {
            if(N_B_integrated[i] > 0.99 * N_B)
                {i_B = i; std::cout << "found i_B at: " << i_B << std::endl;}
        }
        if(i_F < 0) {
            if(N_F_integrated[i] > 0.99 * N_F)
                {i_F = i; std::cout << "found i_F at: " << i_F << std::endl;}
        }
        if(i_B > 0 && i_F > 0)
            break;
    }
    // obtain radius from corresponding index
    double R_B = r[i_B], R_F = r[i_F];

    // compute the fermionic radius R_f using the definition where P(R_f)==0:
    // iterate the Pressure-array until we find the first point where the pressure is zero:
    double R_F_0 = 0.0; // fermionic radius where pressure is zero
    int i_F_0 = 0; // index of fermionic radius wher repressure is zero

    // find the first point where the pressure is (approx) zero and take this as the fermionic radius
    for (unsigned i = 1; i < results.size(); ++i) {
        if (results[i].second[4] < 1e-15) {
            i_F_0 = i;
            break;
        }
    }
    R_F_0 = r[i_F_0];

    std::cout << "M_T = " << M_T << ", N_B = " << N_B << ", R_B = " << R_B << ", N_F = " << N_F << ", R_F = " << R_F << ", R_F_0 = " << R_F_0 << ", N_B/N_F = " << N_B / N_F << std::endl;
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
                << fbs.omega                 << " "   // omega
                << fbs.mu                    << " "   // mass mu
                << fbs.lambda;      // self-interaction parameter lambda
}

