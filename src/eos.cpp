#include "eos.hpp"

/******************
 * PolytropicEoS *
 ******************/
double PolytropicEoS::get_P_from_rho(const double rho, const double epsilon) {
	return this->kappa*std::pow(rho, this->Gamma);
}

double PolytropicEoS::dP_de(const double rho, const double epsilon) { // TODO: Check
	return this->kappa*this->Gamma*std::pow(rho, this->Gamma-1.);
}
/* This expects as input P and will give rho, epsilon as output through reference
 * according to the polytropic EoS
 * TODO : update restmass density (rho) and specific energy (epsilon) according to polytropic EOS */
void PolytropicEoS::callEOS(double& rho, double& epsilon, const double P) {
    rho = std::pow(P / this->kappa, 1./this->Gamma);
    epsilon = this->kappa*std::pow(rho, this->Gamma - 1.) / (this->Gamma - 1.);
    return;
}

double PolytropicEoS::min_P() {
    return 0.;
}
double PolytropicEoS::min_rho() {
    return 0.;
}


/******************
 *  CausalEoS *
 ******************/
double CausalEoS::get_P_from_rho(const double rho, const double epsilon) {

	return this->P_f + rho*(1. + epsilon) - this->eps_f;   // p = P_f + eps - eps_f
}

double CausalEoS::dP_de(const double rho, const double epsilon) {

	return 1.;   // p = P_f + eps - eps_f
}

/* This expects as input P and will give rho, epsilon as output through reference
 * according to the causal EoS
 * TODO : update restmass density (rho) and specific energy (epsilon) according to causal EOS */
void CausalEoS::callEOS(double& rho, double& epsilon, const double P) {
        rho = P - this->P_f + this->eps_f;
        epsilon = 0.;
}

double CausalEoS::min_P() {
    return 0.;
}
double CausalEoS::min_rho() {
    return 0.;
}

/******************
 *  EoStable *
 ******************/


/* This function expects a filename to a table.
 * The table has to contain the rest mass density rho, pressure P, energy density e
 * The indices of these values is given by the std::map such that
 *       indices["rho"] corresponds to the column of rho [1/fm^3]
 *       indices["e"] corresponds to the column of energy density [MeV/fm^3]
 *       indices["P"] corresponds to the column of pressure [MeV/fm^3]
 *
 * Columns are counted starting from 0
 * Lines containing # are ignored
 */
bool EoStable::load_from_file(const std::string filename, std::map<std::string, int> indices) {

    const double MeV_fm3_to_codeunits =	2.886376934e-6; // unit conversion from Mev/fm^3 to code units M_s*c^2/ (G*M_s/c^2)^3
    const double neutron_mass = 939.565379;	 // nuclear density in MeV

    std::ifstream infile;
    infile.open(filename);
    if(!infile.is_open())
        return false;

    rho_table.clear(); P_table.clear(); e_table.clear();
    int max_index = 0.;
    max_index = std::max( std::max(indices.at("rho"), indices.at("P")), indices.at("e"));

    std::cout << "loading EoS table " << filename << std::endl; // TODO: Check error that arises when this is commented out xD
    std::string line;
    while (std::getline(infile, line)) {
        double temp;
        std::stringstream ss(line);

		// ignore lines with '#' in it (to make comment in the table)
		if (line.find('#') != std::string::npos)
    		continue;

        for(int i = 0; i < max_index+1; i++) {
            ss >> temp;
            if( i == indices.at("rho"))
                rho_table.push_back(temp * MeV_fm3_to_codeunits * neutron_mass);
            if( i ==  indices.at("P"))
                P_table.push_back(temp * MeV_fm3_to_codeunits);
            if( i == indices.at("e"))
                e_table.push_back(temp * MeV_fm3_to_codeunits);
        }
    }
    infile.close();
    return true;
}

/* This function expects a filename to a table with the following structure
 *
 * | restmass density rho [1/fm^3] |  (skipped)  | energy density[MeV/fm^3]  |  pressure[MeV/fm^3] |
 *
 * and calls load_from_file with these indices
 */
bool EoStable::load_from_file(const std::string filename) {
    std::map<std::string, int> indices;
    indices["rho"] = 0;
    indices["e"] = 2;
    indices["P"] = 3;
    return load_from_file(filename, indices);
}

/* This expects as input P and will give rho, epsilon as output through reference
 * according to the tabulated EoS
 * If we are below or above the pressures in Pres, 0. is returned
 *  otherwise simple linear interpolation is used to obtain rho, epsilon */
void EoStable::callEOS(double& rho, double& epsilon, const double P) {

	// use P to scan the EOS table (iterate from the top first because we start at high pressures):
	unsigned int table_len = P_table.size();
	double e_tot_tmp = 0.0;


	if (P < P_table[0]) { // we are outside the validity range of the table. Return zero
		rho = 0.;
		epsilon = 0.;
		return;
	}

	// search the table from the smaller values on first
	for (unsigned int i = 1; i<table_len; i++) {
		if (P_table[i] > P) {
			// the correct value is between the i-1th index and the ith index:
			// interpolate linearily between them:
			rho = rho_table[i-1] + (rho_table[i] - rho_table[i-1]) / (P_table[i] - P_table[i-1]) * (P - P_table[i-1]);
			e_tot_tmp = e_table[i-1] + (e_table[i] - e_table[i-1]) / (P_table[i] - P_table[i-1]) * (P - P_table[i-1]);
			epsilon = e_tot_tmp/rho - 1.0;	// re-arrange to get epsilon. e=rho*(1+epsilon)
			return;
		}
	}
	// return 0. outside the validity
    rho = 0.;
    epsilon = 0.;
}

/* This expects as input rho and returns P (epsilon is ignored)
 * according to the tabulated EoS
 * If we are below or above the densities in rho, 0. is returned
 *  otherwise simple linear interpolation is used to obtain P from rho */
double EoStable::get_P_from_rho(const double rho, const double epsilon) {
    unsigned int table_len = rho_table.size();

    // if we are below the table return 0.
    if (rho < rho_table[0])
        return 0.;

    // search the table for the correct value
    for (unsigned int i = 1; i<table_len; i++) {
        if (rho_table[i] > rho) {
            // the correct value is between the ith index and the i+1th index:
            // interpolate linearily between them:
            return P_table[i-1] + (P_table[i] - P_table[i-1]) / (rho_table[i] - rho_table[i-1]) * (rho - rho_table[i-1]);
        }
    }
    // if no value was found return 0.
    return 0.;
}

/* This expects as input rho, epsilon and returns dP/de
 * according to the tabulated EoS
 * If we are below or above the energies in e, 0. is returned
 *  otherwise simple linear interpolation is used to obtain P from e */
double EoStable::dP_de(const double rho, const double epsilon) {
    // search the table for the correct value:
    unsigned int table_len = e_table.size();
    const double e = rho*(1. + epsilon);

    if(e < e_table[1])
        return 0.;

    for (unsigned int i = 2; i<table_len; i++) {
        // scan for the first matching e
        if (e_table[i] > e) {
            // estimate derivative at point i-1 and i
            double dP1 = (P_table[i]-P_table[i-1])/(e_table[i]-e_table[i-1])/2. + (P_table[i-1] - P_table[i-2])/(e_table[i-1]-e_table[i-2])/2.;
            double dP2 = (P_table[i+1]-P_table[i])/(e_table[i+1]-e_table[i])/2. + (P_table[i] - P_table[i-1])/(e_table[i]-e_table[i-1])/2.;
            return dP1 + (dP2-dP1)/(e_table[i]- e_table[i-1]) * (e - e_table[i-1]);
        }
    }
    return 0.;
}

double EoStable::min_P() {
    return this->P_table.at(0);
}
double EoStable::min_rho() {
    return this->rho_table.at(0);
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
