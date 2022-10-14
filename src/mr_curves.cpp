
#include "mr_curves.hpp"


void write_MRphi_curve(const std::vector<FermionBosonStar>& MRphi_curve, std::string filename) {

    std::ofstream img;
	img.open(filename);
    std::vector<std::string> labels = MRphi_curve.at(0).labels();

	if(img.is_open()) {
        // print the labels in line 1:
        img << "# "; // hashtag so that python recognizes it as a commented line
        for(unsigned i = 0; i < labels.size(); i++)
            img << labels[i] << "\t ";
		img << std::endl;

        // print all the data:
        for(auto it = MRphi_curve.begin(); it != MRphi_curve.end(); ++it)
            img << std::scientific << std::setprecision(10) << *it << std::endl;
	}
	img.close();
}

void write_MRphik2_curve(const std::vector<FermionBosonStar>& MRphi_curve,  std::vector<FermionBosonStarTLN>& MRphik2_curve, std::string filename) {

    MRphik2_curve.clear();  MRphik2_curve.reserve(MRphi_curve.size());

    double phi_1_0, phi_1_1;

    time_point start3{clock_type::now()};
    for(auto it = MRphi_curve.begin(); it != MRphi_curve.end(); ++it) {
        FermionBosonStarTLN fbstln(*it);
        phi_1_0 = 1e-3 * it->phi_0;
        phi_1_1 = 1e6 * it->phi_0;
        fbstln.bisection_phi_1(phi_1_0, phi_1_1);
        fbstln.evaluate_model();
        MRphik2_curve.push_back(fbstln);
    }
    time_point end3{clock_type::now()};
    std::cout << "evaluation of "<< MRphik2_curve.size() <<" TLN stars took " << std::chrono::duration_cast<second_type>(end3-start3).count() << "s" << std::endl;
    std::cout << "average time per evaluation: " << (std::chrono::duration_cast<second_type>(end3-start3).count()/(MRphi_curve.size())) << "s" << std::endl;

    if(filename.empty())
        return;

    std::ofstream img;
	img.open(filename);
    std::vector<std::string> labels = MRphik2_curve.at(0).labels();

	if(img.is_open()) {
        // print the labels in line 1:
        img << "# "; // hashtag so that python recognizes it as a commented line
        for(unsigned i = 0; i < labels.size(); i++)
            img << labels[i] << "\t ";
		img << std::endl;

        // print all the data:
        for(auto it = MRphik2_curve.begin(); it != MRphik2_curve.end(); ++it)
            img << std::scientific << std::setprecision(10) << *it << std::endl;
	}
	img.close();

}


// compute curves of constant rho_c and pyh_c:
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

