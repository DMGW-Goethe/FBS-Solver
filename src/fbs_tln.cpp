#include "fbs_tln.hpp"


/* This function gives the initial conditions for the FBS integration
 * with a_0 = 1,  alpha_0 = 1,  phi_0 = this->phi_0, Psi = 0, P_0 = EOS(rho_0)
 *  H_0 = H_0 * r_0^2,  dH_0 = 2 H_0 r , phi_1_0 = phi_0 r_0^3,  dphi_1_0 = phi_0 3 r_0^2
 * at the position r_0=r_init */
vector FermionBosonStarTLN::get_initial_conditions(double r_init) const {
    r_init = (r_init < 0. ? this->r_init : r_init);
    return vector( {1.0, 1.0, this->phi_0, 0., this->rho_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(this->rho_0, 0.) : 0.,
                                                             this->H_0*r_init*r_init, 2.*this->H_0*r_init, this->phi_1_0*pow(r_init,3), 3.*this->phi_1_0*pow(r_init,2)});
}

/* This function gives the system of ODEs for the FBS + TLN star
 *  for the variables a, alpha, phi, Psi, P, H, dH, phi_1, dphi_1
 *
 *  This function is called by the integrator during the integration
 * */
vector FermionBosonStarTLN::dy_dr(const double r, const vector& vars) const {
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
    bool phi_converged = this->phi_0 <= 0.;

    if (this->phi_0 > 0.) {
        if(this->R_B_0 >  0.) // we have artifically set phi to 0 at some point which makes our lifes much easier
            phi_converged = true;
    }
    // std::cout << "calculate_star_parameters with phi_converged = " << phi_converged << std::endl;

    /* The quantity to compute is y = r H' / H
     * at the point where both components have converged */
    auto M_func = [&results](int index) { return results[index].first / 2. * (1. - 1./pow(results[index].second[0], 2)); };
    auto y_func = [&results](int index) { return results[index].first * results[index].second[6]/ results[index].second[5]; };
    auto dy_func = [&results, &y_func] (int i) { return (y_func(i+1) - y_func(i))/(results[i+1].first - results[i].first)/2. + (y_func(i) - y_func(i-1))/(results[i].first - results[i-1].first)/2.;  };
    double y = 0., R_ext = 0., M_ext = 0.;


    if(phi_converged) {
        int index_ext = 0;
        R_ext = std::max(this->R_B_0, this->R_F);
        while(results[index_ext].first < R_ext && index_ext < step_number-2)
            index_ext++;

        y = y_func(index_ext);
        M_ext = M_func(index_ext);
        R_ext = results[index_ext].first;
    }
    else {
        if(this->R_F > 100.*this->R_B) { // in this case, the fermionic part dominates the bosonic part, so we can read out the TLN where the fermionic part vanishes
            int index_R_F = 0;
            while ( results[index_R_F].first < this->R_F && index_R_F < step_number-1)
                index_R_F++;
            // approximate y, M_ext at R_F
            y = y_func(index_R_F-1) +  (y_func(index_R_F) - y_func(index_R_F-1)) / (results[index_R_F].first - results[index_R_F-1].first) * (this->R_F - results[index_R_F-1].first);
            M_ext = M_func(index_R_F-1) +  (M_func(index_R_F) - M_func(index_R_F-1)) / (results[index_R_F].first - results[index_R_F-1].first) * (this->R_F - results[index_R_F-1].first);
            R_ext = R_F;
            //std::cout << "R_F > R_B:  at R_F=" << R_F << " y = " << y << std::endl;
        }
        else { // implement procedure from Sennet 2017 (https://arxiv.org/pdf/1704.08651.pdf)
            // to find the starting point see where y actually has a minimum
            int index_bs_radius = 1;
            while(results[index_bs_radius].first < this->R_B/1e3  && index_bs_radius < step_number-2)
                index_bs_radius++;
            int index_y_min = index_bs_radius;
            while(y_func(index_y_min) < y_func(index_y_min-1) && index_y_min < step_number-2)
                index_y_min++;

            // now look for the local maxima&saddle points of y going from left to right (low r to higher r)
            std::vector<int> indices_maxima;
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
        }
    }
    // std::cout << "R_F=" << R_F <<", R_B=" <<  R_B << ", R_ext = " << R_ext << ", y = " << y << std::endl;

    // now that we found y, calculate k2
    double C = M_ext / R_ext; // the compactness at the extraction point

    /* tidal deformability as taken from https://arxiv.org/pdf/0711.2420.pdf */
    double lambda_tidal = 16./15.*pow(this->M_T, 5)* /* pow(M_ext, 5) */ pow(1.-2.*C, 2)* (2. + 2.*C*(y-1.) - y)
                    / (2.*C*(6. - 3.*y + 3.*C*(5.*y-8.))
                        + 4.*pow(C,3)*(13. - 11.*y + C*(3.*y-2.) + 2.*C*C*(1. + y))
                        + 3.* pow(1. - 2.*C, 2) *(2. - y + 2.*C*(y-1))*  log(1.-2.*C)   );
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
    const bool force_phi_to_0 = true;

    std::vector<integrator::Event> events;
    results.clear();

    std::vector<int> additional_zero_indices = {7,8};
    int res = this->integrate_and_avoid_phi_divergence(results, events,  intOpts, force_phi_to_0, additional_zero_indices);

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
const integrator::Event FermionBosonStarTLN::dphi_1_diverging = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return (std::abs(y[8]) > 1e6); }, true, "dphi_1_diverging");

/* This event triggers when phi_1 is negative, used for zero counting */
const integrator::Event FermionBosonStarTLN::phi_1_negative = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[7] < 0.; }, false, "phi_1_negative");
/* This event triggers when phi_1 is positive  */
const integrator::Event FermionBosonStarTLN::phi_1_positive = integrator::Event([](const double r, const double dr, const vector& y, const vector& dy, const void*params) { return y[7] > 0.; }, false, "phi_1_positive");


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
    const bool force_phi_to_0 = true;

    // variables regarding the integration
    integrator::IntegrationOptions intOpts;
    intOpts.verbose = verbose - 1;
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
    int res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_0, events, intOpts, force_phi_to_0);
    n_roots_0 = events[0].steps.size() + events[1].steps.size() - 1;    // number of roots is number of - to + crossings plus + to - crossings

    // set the upper phi_1 and integrate the ODEs:
    this->phi_1_0 = phi_1_0_r;
    res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_1, events, intOpts, force_phi_to_0);
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
        res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_mid, events, intOpts, force_phi_to_0);
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
    res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_0, events, intOpts, force_phi_to_0);
    n_inft_0 = results_0[results_0.size()-1].second[index_phi_1] > 0.;    // save if sign(Phi_1(inf)) is positive or negative

    this->phi_1_0 = phi_1_0_r;
    res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_1, events, intOpts, force_phi_to_0);
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
        res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_mid, events, intOpts, force_phi_to_0);
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
    res = FermionBosonStar::integrate_and_avoid_phi_divergence(results_0, events, intOpts, force_phi_to_0);

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
