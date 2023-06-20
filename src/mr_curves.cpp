#include "mr_curves.hpp"

/*
Function to write the FBS data into a txt file.
Templated functions must be defined in the hpp file, therefore the following lines are just a placeholder for the sake of readability.
If you want to see the function definition, refer to mr_curves.hpp.
template <typename T>
void write_MRphi_curve(const std::vector<T>& MRphi_curve, std::string filename);
*/

/*
function that calculates all FBS/FPS in a rho-phi grid. Adapted to be able to use both FBS and FPS.
Templated functions must be defined in the hpp file, therefore the following lines are just a placeholder for the sake of readability.
If you want to see the function definition, refer to mr_curves.hpp.
template <typename U>
void calc_rhophi_curves(std::vector<U>& MRphi_curve, int verbose, int mode = 0);
*/

/*
compute curves of constant rho_c and Nb/(Nf+Nb)-ratio:
Templated functions must be defined in the hpp file, therefore the following lines are just a placeholder for the sake of readability.
If you want to see the function definition, refer to mr_curves.hpp.
template <typename V>
void calc_NbNf_curves(NUMERIC mu, NUMERIC lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<NUMERIC>& rho_c_grid, const std::vector<NUMERIC>& NbNf_grid, std::vector<V>& MRphi_curve, int verbose) {
*/

// compute curves of constant rho_c and phi_c.
void calc_rhophi_curves(NUMERIC mu, NUMERIC lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<NUMERIC>& rho_c_grid, const std::vector<NUMERIC>& phi_c_grid, std::vector<FermionBosonStar>& MRphi_curve, int verbose) {

    FermionBosonStar fbs_model(EOS, mu, lambda, 0._num);    // create model for star
    MRphi_curve.clear();
    MRphi_curve.reserve(rho_c_grid.size()*phi_c_grid.size());

    // set initial conditions for every star in the list:
    for(unsigned int j = 0; j < phi_c_grid.size(); j++) {
        for(unsigned int i = 0; i < rho_c_grid.size(); i++) {
            FermionBosonStar fbs(fbs_model);
            fbs.rho_0 = rho_c_grid[i]; fbs.phi_0 = phi_c_grid[j];
            MRphi_curve.push_back(fbs);
        }
    }

    calc_rhophi_curves(MRphi_curve, verbose);
}


// compute the tidal love number for curves of constant rho_c and phi_c:
// the calculation of the unperturbed solution must be performed before, and only then this function can be called because it uses the equilibrium results from calc_rhophi_curves()
void calc_MRphik2_curve(const std::vector<FermionBosonStar>& MRphi_curve,  std::vector<FermionBosonStarTLN>& MRphik2_curve, int verbose) {

	MRphik2_curve.clear();  MRphik2_curve.reserve(MRphi_curve.size());

    // set initial conditions for every star in the list:
    for(auto it = MRphi_curve.begin(); it != MRphi_curve.end(); ++it) {
        FermionBosonStarTLN fbs(*it);	// set the initial conditions for the FBS with tidal love number
        MRphik2_curve.push_back(fbs);
    }

    NUMERIC phi_1_0, phi_1_1; // upper and lower bound for the bisection of the perturbed Phi-field

    unsigned int done = 0;

    time_point start3{clock_type::now()};
	#pragma omp parallel for schedule(dynamic, 10)
    for(unsigned int i = 0; i < MRphi_curve.size(); i++) {
        //FermionBosonStarTLN fbstln(*it);
        phi_1_0 = 1e-3_num * MRphi_curve[i].phi_0;
        phi_1_1 = 1e6_num * MRphi_curve[i].phi_0;
        int bisection_success = MRphik2_curve[i].bisection_phi_1(phi_1_0, phi_1_1);
        if (bisection_success == 0)
            MRphik2_curve[i].evaluate_model();
        //MRphik2_curve[i].push_back(fbstln);

        #pragma omp atomic
        done++;

        if(verbose >1)
            std::cout << "Progress: "<< float(done) / MRphi_curve.size() * 100.0 << "%" << std::endl;
    }

    time_point end3{clock_type::now()};
    if(verbose > 0) {
        std::cout << "evaluation of "<< MRphik2_curve.size() <<" TLN stars took " << std::chrono::duration_cast<second_type>(end3-start3).count() << "s" << std::endl;
        std::cout << "average time per evaluation: " << (std::chrono::duration_cast<second_type>(end3-start3).count()/(MRphi_curve.size())) << "s" << std::endl;
    }
}

void calc_twofluidFBS_curves(std::shared_ptr<EquationOfState> EOS1, std::shared_ptr<EquationOfState> EOS2, const std::vector<NUMERIC>& rho1_c_grid, const std::vector<NUMERIC>& rho2_c_grid, std::vector<TwoFluidFBS>& MRphi_curve, NUMERIC mu, NUMERIC lambda) {
	
	TwoFluidFBS fbs_model(EOS1, EOS2, mu, lambda);    // create model for star
	
    MRphi_curve.clear();
    MRphi_curve.reserve(rho1_c_grid.size()*rho2_c_grid.size());

    // set initial conditions for every star in the list:
    for(unsigned int j = 0; j < rho2_c_grid.size(); j++) {
        for(unsigned int i = 0; i < rho1_c_grid.size(); i++) {
            TwoFluidFBS fbs(fbs_model);
            fbs.rho1_0 = rho1_c_grid[i]; fbs.rho2_0 = rho2_c_grid[j];
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


// compute curves of constant rho_c and E_c for the Fermion-Proca Star:
void calc_rhophi_curves_FPS(NUMERIC mu, NUMERIC lambda, std::shared_ptr<EquationOfState> EOS, const std::vector<NUMERIC>& rho_c_grid, const std::vector<NUMERIC>& E_c_grid, std::vector<FermionProcaStar>& MRphi_curve, int verbose, int mode) {

    FermionProcaStar fps_model(EOS, mu, lambda);    // create model for star
    MRphi_curve.clear();
    MRphi_curve.reserve(rho_c_grid.size()*E_c_grid.size());

    // set initial conditions for every star in the list:
    for(unsigned int j = 0; j < E_c_grid.size(); j++) {
        for(unsigned int i = 0; i < rho_c_grid.size(); i++) {
            FermionProcaStar fps(fps_model);
            fps.rho_0 = rho_c_grid[i]; fps.E_0 = E_c_grid[j]; fps.phi_0 = fps.E_0;
            MRphi_curve.push_back(fps);
        }
    }

    calc_rhophi_curves(MRphi_curve, verbose, mode);
}