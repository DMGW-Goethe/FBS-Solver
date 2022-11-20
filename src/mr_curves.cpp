
#include "mr_curves.hpp"

/*
Function to write the FBS data into a txt file.
Templated functions must be defined in the hpp file, therefore the following lines are just a placeholder for the sake of readability.
If you want to see the function definition, refer to mr_curves.hpp.
template <typename T>
void write_MRphi_curve(const std::vector<T>& MRphi_curve, std::string filename);
*/

// compute curves of constant rho_c and phi_c:
void calc_rhophi_curves(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& phi_c_grid, std::vector<FermionBosonStar>& MRphi_curve) {

    FermionBosonStar fbs_model(EOS, mu, lambda, 0.);    // create model for star
    MRphi_curve.clear();
    MRphi_curve.reserve(rho_c_grid.size()*phi_c_grid.size());

    // set initial conditions for every star in the list:
    for(unsigned int j = 0; j < phi_c_grid.size(); j++) {
        for(unsigned int i = 0; i < rho_c_grid.size(); i++) {
            FermionBosonStar fbs(fbs_model);
            fbs.set_initial_conditions(rho_c_grid[i], phi_c_grid[j]);
            MRphi_curve.push_back(fbs);
        }
    }

    // calculate curve
    const double omega_0 = 1., omega_1 = 10.;  // upper and lower bound for omega in the bisection search

    time_point start3{clock_type::now()};
    #pragma omp parallel for
    for(unsigned int i = 0; i < MRphi_curve.size(); i++) {
        int bisection_success = MRphi_curve[i].bisection(omega_0, omega_1);  // compute bisection
		if (bisection_success == -1)
            std::cout << "Bisection failed with omega_0=" << omega_0 << ", omega_1=" << omega_1 << " for " << MRphi_curve[i] << std::endl;
        MRphi_curve[i].evaluate_model();   // evaluate the model but do not save the intermediate data into txt file
    }
    time_point end3{clock_type::now()};
    std::cout << "evaluation of "<< MRphi_curve.size() <<" stars took " << std::chrono::duration_cast<second_type>(end3-start3).count() << "s" << std::endl;
    std::cout << "average time per evaluation: " << (std::chrono::duration_cast<second_type>(end3-start3).count()/(MRphi_curve.size())) << "s" << std::endl;

}


// compute curves of constant rho_c and Nb/NF-ratio:
void calc_NbNf_curves(double mu, double lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<double>& rho_c_grid, const std::vector<double>& NbNf_grid, std::vector<FermionBosonStar>& MRphi_curve) {
    MRphi_curve.clear();
    MRphi_curve.reserve(rho_c_grid.size()*NbNf_grid.size());

    // set initial conditions for every star in the list:
    for (unsigned j = 0; j < NbNf_grid.size(); j++) {
        for(unsigned i = 0; i < rho_c_grid.size(); i++) {
            FermionBosonStar fbs(EOS, mu, lambda, 0.);  // create a star
            fbs.set_initial_conditions(rho_c_grid[i], 1e-10);
            MRphi_curve.push_back(fbs);
        }
    }

    // compute the MR-diagrams:
    double omega_0 = 1., omega_1 = 10.;

    #pragma omp parallel for
    for(unsigned i = 0; i < NbNf_grid.size(); i++) {
        for(unsigned j = 0; j < rho_c_grid.size() ; j++) {
            int index = i*rho_c_grid.size() + j;
            MRphi_curve[index].shooting_NbNf_ratio(NbNf_grid[i], 1e-4, omega_0, omega_1);  // compute star with set NbNf ratio
            //MRphi_curve[index].evaluate_model();   // evaluate the model but do not save the intermediate data into txt file

        }
    }

    std::cout << "end loop" << std::endl;
}


// compute the tidal love number for curves of constant rho_c and phi_c:
// the calculation of the unperturbed solution must be performed before, and only then this function can be called because it uses the equilibrium results from calc_rhophi_curves()
void calc_MRphik2_curve(const std::vector<FermionBosonStar>& MRphi_curve,  std::vector<FermionBosonStarTLN>& MRphik2_curve) {

	MRphik2_curve.clear();  MRphik2_curve.reserve(MRphi_curve.size());

    // set initial conditions for every star in the list:
    for(auto it = MRphi_curve.begin(); it != MRphi_curve.end(); ++it) {
        FermionBosonStarTLN fbs(*it);	// set the initial conditions for the FBS with tidal love number
        MRphik2_curve.push_back(fbs);
    }

    double phi_1_0, phi_1_1; // upper and lower bound for the bisection of the perturbed Phi-field

    time_point start3{clock_type::now()};
	#pragma omp parallel for
    for(unsigned int i = 0; i < MRphi_curve.size(); i++) {
        //FermionBosonStarTLN fbstln(*it);
        phi_1_0 = 1e-3 * MRphi_curve[i].phi_0;
        phi_1_1 = 1e6 * MRphi_curve[i].phi_0;
        MRphik2_curve[i].bisection_phi_1(phi_1_0, phi_1_1);
        MRphik2_curve[i].evaluate_model();
        //MRphik2_curve[i].push_back(fbstln);
    }

    time_point end3{clock_type::now()};
    std::cout << "evaluation of "<< MRphik2_curve.size() <<" TLN stars took " << std::chrono::duration_cast<second_type>(end3-start3).count() << "s" << std::endl;
    std::cout << "average time per evaluation: " << (std::chrono::duration_cast<second_type>(end3-start3).count()/(MRphi_curve.size())) << "s" << std::endl;
}

//
void calc_twofluid_curves(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2, const std::vector<double>& rho1_c_grid, const std::vector<double>& rho2_c_grid, std::vector<TwoFluidFBS>& MRphi_curve) {

	TwoFluidFBS fbs_model(EOS1, EOS2);    // create model for star
    MRphi_curve.clear();
    MRphi_curve.reserve(rho1_c_grid.size()*rho2_c_grid.size());

    // set initial conditions for every star in the list:
    for(unsigned int j = 0; j < rho2_c_grid.size(); j++) {
        for(unsigned int i = 0; i < rho1_c_grid.size(); i++) {
            TwoFluidFBS fbs(fbs_model);
            fbs.set_initial_conditions(rho1_c_grid[i], rho2_c_grid[j]);
            MRphi_curve.push_back(fbs);
        }
    }

	time_point start3{clock_type::now()};
	// integrate all the stars in parallel:
    #pragma omp parallel for
    for(unsigned int i = 0; i < MRphi_curve.size(); i++) {
        MRphi_curve[i].evaluate_model();   // evaluate the model but do not save the intermediate data into txt file
    }
    time_point end3{clock_type::now()};
    std::cout << "evaluation of "<< MRphi_curve.size() <<" stars took " << std::chrono::duration_cast<second_type>(end3-start3).count() << "s" << std::endl;
    std::cout << "average time per evaluation: " << (std::chrono::duration_cast<second_type>(end3-start3).count()/(MRphi_curve.size())) << "s" << std::endl;

}