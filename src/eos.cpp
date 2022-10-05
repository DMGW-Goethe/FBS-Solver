#include "eos.hpp"

/* PolytropicEoS */
double PolytropicEoS::get_P_from_rho(const double rho_in, const double epsilon) {

	return this->kappa*std::pow(rho_in, this->Gamma);	// p = k*rho^gamma
}

double PolytropicEoS::dP_drho(const double rho_in, const double epsilon) {

	return this->kappa*this->Gamma*std::pow(rho_in, this->Gamma-1.);
}

void PolytropicEoS::callEOS(double& myrho, double& epsilon, const double P) {
	// update restmass density (rho) and specific energy (epsilon) according to polytropic EOS
    	myrho = std::pow(P / this->kappa, 1./this->Gamma);
    	epsilon = this->kappa*std::pow(myrho, this->Gamma - 1.) / (this->Gamma - 1.);
		return;
}

double PolytropicEoS::min_P() {
    return 0.;
}
double PolytropicEoS::min_rho() {
    return 0.;
}


/* CausalEoS */
double CausalEoS::get_P_from_rho(const double rho_in, const double epsilon) {

	return this->P_f + rho_in*(1. + epsilon) - this->eps_f;   // p = P_f + eps - eps_f
}

double CausalEoS::dP_drho(const double rho_in, const double epsilon) {

	return (1. + epsilon);   // p = P_f + eps - eps_f
}

void CausalEoS::callEOS(double& myrho, double& epsilon, const double P) {
	// update restmass density (rho) and specific energy (epsilon) according to causal EOS
        myrho = P - this->P_f + this->eps_f;
        epsilon = 0.;
}

double CausalEoS::min_P() {
    return 0.;
}
double CausalEoS::min_rho() {
    return 0.;
}

/* EoStable */
EoStable::EoStable(const std::string filename) {
    // load in the tabulated EOS data and convert it to code units directly

    eos_table_name = filename;
    const double MeV_fm3_to_codeunits =	2.886376934e-6; // unit conversion from Mev/fm^3 to code units M_s*c^2/ (G*M_s/c^2)^3
    const double neutron_mass = 939.565379;	 // nuclear density in MeV

    std::ifstream infile;	// declare input file stream
    infile.open(eos_table_name);

    if(!infile.is_open())
        throw std::runtime_error("File not open");
    std::string line;

    while (std::getline(infile, line)) {

        double temp;
        std::stringstream ss(line);

		// ignore lines with '#' in it (to make comment in the table)
		if (line.find('#') != std::string::npos)
    		continue; // '#' found, therefore we skip this line

		// ---------------------------------------
		// read in the EOS format from the CompOSE code. Format is:
		// n_b[1/fm^3]  Y_e[dimensionless]  energy density[MeV/fm^3]  p[MeV/fm^3]

        ss >> temp; // column 1 -> restmass density
        rho.push_back(temp*MeV_fm3_to_codeunits*neutron_mass);	// write the 1nd column -> restmass density and convert to code units

        ss >> temp; // column 2 -> electron (lepton) fraction (we don't need it, skip this column)
		ss >> temp; // column 3 -> energy density
        e_tot.push_back(temp*MeV_fm3_to_codeunits);	// write the 3nd column -> energy density e=rho*(1+epsilon) and convert to code units

        ss >> temp; // column 4 -> pressure
        Pres.push_back(temp*MeV_fm3_to_codeunits);	// write the 4nd column -> pressure and convert to code units

    }
    infile.close();	// close file after reading from it

}

// call EOS to obtain rho and e from input P:
void EoStable::callEOS(double& myrho, double& epsilon, const double P) {

	// use P to scan the EOS table (iterate from the top first because we start at high pressures):
	unsigned int table_len = Pres.size();
	double e_tot_tmp = 0.0;


	if (P < Pres[0]) { // we are close to vacuum and don't need to search the table. Interpolate with zero
		myrho = 0.0 + (rho[0] / Pres[0]) * P;
		e_tot_tmp = 0.0 + (e_tot[0] / Pres[0]) * P;
		epsilon = e_tot_tmp/myrho - 1.0;	// re-arrange to get epsilon. e=rho*(1+epsilon)
		return;
	}

	// search the table from the smaller values on first
	for (unsigned int i = 1; i<table_len; i++) {
		if (Pres[i] > P) {
			// the correct value is between the i-1th index and the ith index:
			// interpolate linearily between them:
			myrho = rho[i-1] + (rho[i] - rho[i-1]) / (Pres[i] - Pres[i-1]) * (P - Pres[i-1]);
			e_tot_tmp = e_tot[i-1] + (e_tot[i] - e_tot[i-1]) / (Pres[i] - Pres[i-1]) * (P - Pres[i-1]);
			epsilon = e_tot_tmp/myrho - 1.0;	// re-arrange to get epsilon. e=rho*(1+epsilon)
			return;
		}
	}
	// extrapolate linearly
	myrho = rho[table_len-2] + (rho[table_len-1] - rho[table_len-2]) / (Pres[table_len-1] - Pres[table_len-2]) * (P - Pres[table_len-2]);
	e_tot_tmp = e_tot[table_len-2] + (e_tot[table_len-1] - e_tot[table_len-2]) / (Pres[table_len-1] - Pres[table_len-2]) * (P - Pres[table_len-2]);
	epsilon = e_tot_tmp/myrho - 1.0;
}

// obtain P from a rho input
double EoStable::get_P_from_rho(const double rho_in, const double epsilon) {
	unsigned int table_len = rho.size();

    // if we are below the table interpolate between (0.,0.) and (rho[0], P[0])
    if (rho_in <= rho[0]) {
        return 0. + (Pres[0]-0.) / (rho[0]-0.) * (rho_in-0.);
    }

	// search the table for the correct value
	for (unsigned int i = 1; i<table_len; i++) {
		if (rho[i] > rho_in) {
			// the correct value is between the ith index and the i+1th index:
			// interpolate linearily between them:
			return Pres[i-1] + (Pres[i] - Pres[i-1]) / (rho[i] - rho[i-1]) * (rho_in - rho[i-1]);
		}
	}
    // if no value was found extrapolate linearly
	return Pres[table_len-2] + (Pres[table_len-1] - Pres[table_len-2]) / (rho[table_len-1] - rho[table_len-2]) * (rho_in - rho[table_len-2]);
}

// obtain dP/drho from a rho input
double EoStable::dP_drho(const double rho_in, const double epsilon) {
	// search the table for the correct value:
	unsigned int table_len = rho.size();

    //assert(rho_in < rho[table_len-2]); // out of range will return 0.
    if(rho_in < rho[1]) {
        double dP1 = (Pres[2]-Pres[1])/(rho[2]-rho[1])/2. + (Pres[1] - Pres[0])/(rho[1]-rho[0])/2.;
        return 0. + dP1 / (rho[1] - 0.) * (rho_in - 0.);
    }

	for (unsigned int i = 2; i<table_len; i++) {
		// scan for the first matching rho
		if (rho[i] > rho_in) {
            // estimate derivative at point i-1 and i
            double dP1 = (Pres[i]-Pres[i-1])/(rho[i]-rho[i-1])/2. + (Pres[i-1] - Pres[i-2])/(rho[i-1]-rho[i-2])/2.;
            double dP2 = (Pres[i+1]-Pres[i])/(rho[i+1]-rho[i])/2. + (Pres[i] - Pres[i-1])/(rho[i]-rho[i-1])/2.;
            return dP1 + (dP2-dP1)/(rho[i]- rho[i-1]) * (rho_in - rho[i-1]);
		}
	}
	return 0.;
}

double EoStable::min_P() {
    return this->Pres.at(0);
}
double EoStable::min_rho() {
    return this->rho.at(0);
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
