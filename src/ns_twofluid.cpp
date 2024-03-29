#include "ns_twofluid.hpp"

using namespace FBS;

/***********************
 * NSTwoFluid *
 ***********************/

// initial conditions for two arbitrary fluids
vector NSTwoFluid::get_initial_conditions(const double r_init) const {

    return vector({0.0, 0.0, 0.0, rho1_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(rho1_0, 0.) : 0.,
                    rho2_0 > this->EOS_fluid2->min_rho() ? this->EOS_fluid2->get_P_from_rho(rho2_0, 0.) : 0., 2.0});

}

vector NSTwoFluid::dy_dr(const double r, const vector &vars) const {

    const double /*nu = vars[0],*/ m1 = vars[1], m2 = vars[2];
    double P1 = vars[3], P2 = vars[4];
    const double Y = vars[5];

    // load the EoS from both fluids:
    EquationOfState &myEOS1 = *(this->EOS);
    EquationOfState &myEOS2 = *(this->EOS_fluid2);

    // call both EoS and compute the wanted values
    double etot1=0., dP1_detot=0., detot_dP1=0.;
    double etot2=0., dP2_detot=0., detot_dP2=0.;
    // first EoS
    if (P1 <= 0. || P1 < myEOS1.min_P())
    {
        P1 = 0.;
        etot1 = 0.;
    }
    else
    {
        etot1 = myEOS1.get_e_from_P(P1);
        dP1_detot = etot1 > myEOS1.min_e() ? myEOS1.dP_de(etot1) : 0.;
        detot_dP1 = dP1_detot > 0. ? 1. / dP1_detot : 0.;
    }
    // second EoS
    if (P2 <= 0. || P2 < myEOS2.min_P())
    {
        P2 = 0.;
        etot2 = 0.;
    }
    else
    {
        etot2 = myEOS2.get_e_from_P(P2);
        dP2_detot = etot2 > myEOS2.min_e() ? myEOS2.dP_de(etot2) : 0.;
        detot_dP2 = dP2_detot > 0. ? 1. / dP2_detot : 0.;
    }

    // compute 'total'-values and helper variables:
    double mtot = m1 + m2;
    double Ptot = P1 + P2;
    double etot_tot = etot1 + etot2;
    double e_lambda = 1. / (1. - 2. * mtot / r);						 // metric component: e^{-lambda(r)} = 1 - 2mtot(r)/r
    double rhoP_param = (etot1 + P1) * detot_dP1 + (etot2 + P2) * detot_dP2; // special parameter used to compute Q

    // compute the ODE:
    double dnu_dr = 2. * (mtot + 4 * M_PI * r * r * r * Ptot) * e_lambda / (r * r);
    double dm1_dr = 4. * M_PI * r * r * etot1;
    double dm2_dr = 4. * M_PI * r * r * etot2;
    double dP1_dr = -0.5 * dnu_dr * (etot1 + P1);
    double dP2_dr = -0.5 * dnu_dr * (etot2 + P2);
    // compute tidal function:
    double Q = 4. * M_PI * e_lambda * (5. * (etot_tot) + 9. * (P1 + P2) + rhoP_param) - 6. * e_lambda / r / r - dnu_dr * dnu_dr; // helper variable Q
    double dY_dr = -Y * Y / r - Y * e_lambda * (1. + 4. * M_PI * r * r * (Ptot - etot_tot)) / r - r * Q;								// tidal perturbation function y(r)

    return vector({dnu_dr, dm1_dr, dm2_dr, dP1_dr, dP2_dr, dY_dr});
}


void NSTwoFluid::calculate_star_parameters(const std::vector<integrator::step> &results, const std::vector<integrator::Event> &events) {

    // compute all values of the star, including the tidal deformability:

    int last_index = results.size() - 1; // last index in the integration

    // next compute the masses of only 1st fluid and 2nd fluid respectively:
    // it is fine to take the last indices of the whole integration because once one fluid has reached zero pressure, the change in mass dm/dr will be zero
    double M_1 = results[last_index].second[1];
    double M_2 = results[last_index].second[2];

    // compute the radii of 1st and 2nd fluid:
    // radius where pressure is approximately zero:
    int i_R1_0 = last_index, i_R2_0 = last_index;
    double R_1_0 = results[last_index].first, R_2_0 = results[last_index].first;
    // there is an edgecase where the pressure does not drop below P_ns_min within the largest allowed integration radius R_MAX (=500).
    // in this case the below iteration will return the initialization value for R_1_0, i_R1_0 etc.
    for (unsigned i = 1; i < results.size(); i++)
    {
        if (results[i].second[3] < P_ns_min || results[i].second[3] < this->EOS->min_P())
        {
            R_1_0 = results[i].first; // 1st fluid
            i_R1_0 = i;
            break;
        }
    }
    for (unsigned i = 1; i < results.size(); i++)
    {
        if (results[i].second[4] < P_ns_min || results[i].second[4] < this->EOS_fluid2->min_P())
        {
            R_2_0 = results[i].first; // 2nd fluid
            i_R2_0 = i;
            break;
        }
    }

    // obtain estimate for the total mass:
    // the total mass M_T is defined at the point M(R) where R ist the largest radius between fluid 1 and fluid 2:
    double M_T = M_1 + M_2; 	// total mass as sum of constituent masses (actually energy-momentum-content) of both fluids

    // radius where 99% of the fluid (1st/2nd respectively) is contained:
    int i_R_1 = 0, i_R_2 = 0;
    for (unsigned int i = 0; i < results.size(); i++)
    {
        if (results[i].second[1] < 0.99 * M_1)
            i_R_1++;
        if (results[i].second[2] < 0.99 * M_2)
            i_R_2++;
    }
    // obtain radius from corresponding index
    double R_1 = results[i_R_1].first, R_2 = results[i_R_2].first;

    // compute the fermion/boson number using the conserved Noether current
    /*  N_B, N_F i.e. N_1 and N_2
     *  We need to integrate the particle number densities to obtain N_B, N_F */
    std::vector<double> r(results.size()), N_B_integrand(results.size()), N_F_integrand(results.size());
    vector v;
    double rho, eps;

    for(unsigned int i = 0; i < results.size(); i++) {
        r[i] = results[i].first;
        v = results[i].second;
        // calc values for 1st fluid: P1 = v[3]
        if (v[3] < P_ns_min || v[3] < this->EOS->min_P()) {rho = 0.;}
        else {this->EOS->callEOS(rho, eps, v[3]);}
        double sqrt_g_rr = std::sqrt( 1./( 1. - 2.*(v[1]+v[2])/r[i] ) );	// calculate the metric component sqrt(g_rr) = sqrt( 1 - 2m(r)/r  )^-1
        N_F_integrand[i] = 4.*M_PI * sqrt_g_rr * rho * r[i] * r[i];   // get fermionic mass (paricle number) for each r

        // do the same for the 2nd fluid: P2 = v[4]
        if (v[4] < P_ns_min || v[4] < this->EOS_fluid2->min_P()) {rho = 0.;}
        else {this->EOS_fluid2->callEOS(rho, eps, v[4]);}
        // sqrt_g_rr is the same for both fluids, so no need to compute it again
        N_B_integrand[i] = 4.*M_PI * sqrt_g_rr * rho * r[i] * r[i]; // get bosonic mass (paricle number) for each r
    }

    // Integrate
    std::vector<double> N_F_integrated, N_B_integrated;
    integrator::cumtrapz(r, N_F_integrand, N_F_integrated);
    integrator::cumtrapz(r, N_B_integrand, N_B_integrated);

    // Compute the conserved Noether charges (i.e. fermion number and boson number)
    double N_F =  N_F_integrated[i_R1_0],
           N_B =  N_B_integrated[i_R2_0];
    // ---------------------------------------------------------------------------------

    // compute the tidal deformability:
    // it is defined at the point where both fluids reach zero density (i.e. the surface of the combined star)
    double maxR = std::max(R_1_0, R_2_0);
    double C = M_T / maxR; // compactness of conbined configuration
    // otain the value of y(r) at maxR:
    int maxRindex = std::max(i_R1_0, i_R2_0);
    double y_R = results[maxRindex].second[5];
    //std::cout << y_R << " " << C << std::endl;

    // compute tidal love number k2:
    // equation taken from: PHYSICAL REVIEW D 105, 123010 (2022)
    // Note: in this model it is possible that some configurations will have a compactness >0.5, meaning that these objects are Black Holes
    // The Buchdal-limit actually allows only stable compact objects with compactness > 4/9 ~ 0.4444. Higher compactness means that the object is not stable and will collapse to a BH
    double k2 = (8. / 5.) * std::pow(C, 5) * std::pow(1. - 2. * C, 2) * (2. - y_R + 2. * C * (y_R - 1.)) /
                (2. * C * (6. - 3. * y_R + 3. * C * (5. * y_R - 8.)) + 4. * std::pow(C, 3) * (13. - 11. * y_R + C * (3. * y_R - 2.) + 2. * C * C * (1. + y_R)) + 3. * std::pow(1. - 2. * C, 2) * (2. - y_R + 2. * C * (y_R - 1.) )* std::log(1. - 2. * C));

    double lambda_tidal = (2. / 3.) * k2 * std::pow(maxR, 5);

    // update all the global star values:
    this->M_T = M_T;
    this->M_1 = M_1;
    this->M_2 = M_2;
    this->R_1 = R_1;
    this->R_1_0 = R_1_0;
    this->R_2 = R_2;
    this->R_2_0 = R_2_0;
    this->C = C;	// compactness
    this->k2 = k2;
    this->lambda_tidal = lambda_tidal;
    this->N_1 = N_F;
    this->N_2 = N_B;
}

void NSTwoFluid::evaluate_model() {

    std::vector<integrator::step> results;
    this->evaluate_model(results);
}

void NSTwoFluid::evaluate_model(std::vector<integrator::step> &results, std::string filename) {

    // define variables used in the integrator and events during integration:
    integrator::IntegrationOptions intOpts;
    intOpts.save_intermediate = true;
    // stop integration if pressure is zero for both fluids:
    std::vector<integrator::Event> events = {NSTwoFluid::all_Pressure_zero};
    results.clear();

    this->integrate(results, events, this->get_initial_conditions(), intOpts); // integrate the star


    // option to save all the radial profiles into a txt file:
    if (!filename.empty())
    {
        plotting::save_integration_data(results, {0, 1, 2, 3, 4, 5}, {"nu", "m1", "m2", "P1", "P2", "y(r)"}, filename);
    }

    this->calculate_star_parameters(results, events);
}

std::ostream& FBS::operator<<(std::ostream &os, const NSTwoFluid &fbs) {

    return os << fbs.M_T << " "					// total gravitational mass
              << fbs.rho1_0 << " "				// central density of 1st fluid
              << fbs.rho2_0 << " "				// central density of 2nd fluid / or central value of scalar field, if effective bosonic EOS is used
              << fbs.R_1 * 1.476625061 << " "	// radius (99% matter included) of 1st fluid
              << fbs.R_1_0 * 1.476625061 << " " // radius where P(r)=0 of 1st fluid
              << fbs.M_1 << " "					// total mass of 1st fluid
              << fbs.N_1 << " "					// total particle number of 1st fluid
              << fbs.R_2 * 1.476625061 << " "	// radius (99% matter included) of 2nd fluid
              << fbs.R_2_0 * 1.476625061 << " " // radius where P(r)=0 of 2nd fluid
              << fbs.M_2 << " "					// total mass of 2nd fluid
              << fbs.N_2 << " "					// total particle number of 2nd fluid
              << fbs.M_2 / fbs.M_1 << " "		// mass ratio M_2 / M_1
              << fbs.N_2 / fbs.N_1 << " "		// particle number ratio N_2 / N_1
              << fbs.C << " "					// Compactness
              << fbs.k2 << " "
              << fbs.lambda_tidal;
}

std::vector<std::string> NSTwoFluid::labels() {
	// labels for the pure two-fluid case:
    return std::vector<std::string>({"M_T", "rho1_0", "rho2_0", "R_1", "R_1_0", "M_1", "N_1", "R_2", "R_2_0", "M_2", "N_2", "M_2/M_1", "N_2/N_1", "C", "k2", "lambda_tidal"});
}

const integrator::Event NSTwoFluid::all_Pressure_zero = integrator::Event([](const double r, const double dr, const vector &y, const vector &dy, const void *params)
                                                                           { return ((y[3] <= 0.1*P_ns_min) && (y[4] <= 0.1*P_ns_min)); },
                                                                           true);
