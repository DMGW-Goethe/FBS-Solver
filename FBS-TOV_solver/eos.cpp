#include "eos.hpp"

/* PolytropicEoS */
double PolytropicEoS::get_P_from_rho(const double rho_in) {

	return this->kappa*std::pow(rho_in, this->Gamma);	// p = k*rho^gamma
}

void PolytropicEoS::callEOS(double& myrho, double& epsilon, const double P) {
    	myrho = std::pow(P / this->kappa, 1./this->Gamma); // update restmass density according to polytropic EOS
    	epsilon = this->kappa*std::pow(myrho, this->Gamma - 1.) / (this->Gamma - 1.);
		return;
}


EoStable::EoStable(const std::string filename) {
    // load in the tabulated EOS data and convert it to code units directly

    eos_table_name = filename;
    const double MeV_fm3_to_codeunits =	2.886376934e-6; // unit conversion from Mev/fm^3 to code units M_s*c^2/ (G*M_s/c^2)^3
    const double neutron_mass = 939.565379;	 // nuclear density in MeV

    std::ifstream infile;	// declare input file stream
    infile.open(eos_table_name);

    assert(infile.is_open());
    std::string line;

    while (std::getline(infile, line)) {

        double temp;
        std::stringstream ss(line);

        ss >> temp; // column 1 -> temperature in MeV (we don't need it)
        ss >> temp; // column 2
        rho.push_back(temp*MeV_fm3_to_codeunits*neutron_mass);	// write the 2nd column -> restmass density and convert to code units

        ss >> temp; // column 3 -> hadronic charge fraction (we don't need it)
        ss >> temp; // column 4
        Pres.push_back(temp*MeV_fm3_to_codeunits);	// write the 4nd column -> pressure and convert to code units

        ss >> temp; // column 5
        e_tot.push_back(temp*MeV_fm3_to_codeunits);	// write the 5nd column -> energy density e=rho*(1+epsilon) and convert to code units
    }
    infile.close();	// close file after reading from it

}

// call EOS to obtain rho and e from input P:
void EoStable::callEOS(double& myrho, double& epsilon, const double P) {

	// use P to scan the EOS table (iterate from the top first because we start at high pressures):
	int table_len = Pres.size();
	double e_tot_tmp = 0.0;

    assert(P >= 0.);
    assert(P < Pres[table_len-1]);

	if (P < Pres[0]) { // we are close to vacuum and don't need to search the table. Interpolate with zero
		myrho = 0.0 + (rho[0] / Pres[0]) * P;
		e_tot_tmp = 0.0 + (e_tot[0] / Pres[0]) * P;
		epsilon = e_tot_tmp/myrho - 1.0;	// re-arrange to get epsilon. e=rho*(1+epsilon)
		return;
	}

	for (int i = table_len-1; i>=0; i--) {
		// scan for the first matching P
		if (Pres[i]< P) {

			// the correct value is between the ith index and the i+1th index:
			// interpolate linearily between them:

			myrho = rho[i] + (rho[i+1] - rho[i]) / (Pres[i+1] - Pres[i]) * (P - Pres[i]);
			e_tot_tmp = e_tot[i] + (e_tot[i+1] - e_tot[i]) / (Pres[i+1] - Pres[i]) * (P - Pres[i]);

			epsilon = e_tot_tmp/myrho - 1.0;	// re-arrange to get epsilon. e=rho*(1+epsilon)

			return;
		}
	}
}

// obtain P from a rho input
double EoStable::get_P_from_rho(const double rho_in) {
	// search the table for the correct value:
	unsigned table_len = rho.size();

    assert(rho_in < rho[table_len-1]);
    assert(rho_in > rho[0]);
	for (int i = table_len-1; i>=0; i--) {
		if (rho[i]< rho_in) {
			// the correct value is between the ith index and the i+1th index:
			// interpolate linearily between them:
			return  Pres[i] + (Pres[i+1] - Pres[i]) / (rho[i+1] - rho[i]) * (rho_in - rho[i]);
		}
	}
	return 0.;
}


/*
#units
uc = 2.99792458*10**(10)	// c_0 in cgs units
uG = 6.67428*10**(-8)		// gravitational constant in cgs units
uMs = 1.9884*10**(33)		// solar mass in cgs units
utime = uG*uMs/uc**3*1000	// time in milliseconds
ulenght = (uG*uMs/uc**2)/100000	// length in km
urho = (uc**6)/(uG**3*uMs**2)	// density in cgs units
normalnuc = 2.705*10**(14)		// normal nuclear density in cgs units
*/
