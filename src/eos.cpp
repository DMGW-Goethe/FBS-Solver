#include "eos.hpp"

/******************
 * PolytropicEoS *
 ******************/
NUMERIC PolytropicEoS::get_P_from_rho(const NUMERIC rho, const NUMERIC epsilon) {
	return this->kappa*std::pow(rho, this->Gamma);
}

NUMERIC PolytropicEoS::get_P_from_e(const NUMERIC etot_in) {
	// todo, maybe have to do root finding?
	throw std::runtime_error("not implemented");
	return 0.0_num;
}

NUMERIC PolytropicEoS::get_e_from_P(const NUMERIC P_in) {
	// etot = rho*(1+epsilon)
	NUMERIC myrho = std::pow(P_in / this->kappa, 1./this->Gamma);
    NUMERIC epsilon = this->kappa*std::pow(myrho, this->Gamma - 1._num) / (this->Gamma - 1._num);
	return ( myrho*(1._num + epsilon) );
}

NUMERIC PolytropicEoS::dP_drho(const NUMERIC rho_in, const NUMERIC epsilon) {

	return this->kappa*this->Gamma*std::pow(rho_in, this->Gamma-1._num);
}

NUMERIC PolytropicEoS::dP_de(const NUMERIC e) {
	throw std::runtime_error("not implemented");
	return 0.0_num;	// TBD
}


/* This expects as input P and will give rho, epsilon as output through reference
 * according to the polytropic EoS
 * TODO : update restmass density (rho) and specific energy (epsilon) according to polytropic EOS */
void PolytropicEoS::callEOS(NUMERIC& rho, NUMERIC& epsilon, const NUMERIC P) {
    rho = std::pow(P / this->kappa, 1._num/this->Gamma);
    epsilon = this->kappa*std::pow(rho, this->Gamma - 1._num) / (this->Gamma - 1._num);
    return;
}

NUMERIC PolytropicEoS::min_P() {
    return 0._num;
}
NUMERIC PolytropicEoS::min_rho() {
    return 0._num;
}
NUMERIC PolytropicEoS::min_e() {
    return 0._num;
}


/******************
 *  CausalEoS *
 ******************/
NUMERIC CausalEoS::get_P_from_rho(const NUMERIC rho, const NUMERIC epsilon) {

	return this->P_f + rho*(1._num + epsilon) - this->eps_f;   // p = P_f + eps - eps_f
}


NUMERIC CausalEoS::get_P_from_e(const NUMERIC etot) {

    return this->P_f + etot - this->eps_f;   // p = P_f + eps - eps_f
}

NUMERIC CausalEoS::get_e_from_P(const NUMERIC P) {

    return P - this->P_f + this->eps_f;   // p = P_f + eps - eps_f
}

NUMERIC CausalEoS::dP_de(const NUMERIC e) {

	return 1.;   // p = P_f + eps - eps_f
}

NUMERIC CausalEoS::dP_drho(const NUMERIC rho_in, const NUMERIC epsilon) {

	return 1._num+epsilon;
}

/* This expects as input P and will give rho, epsilon as output through reference
 * according to the causal EoS
 * TODO : update restmass density (rho) and specific energy (epsilon) according to causal EOS */
void CausalEoS::callEOS(NUMERIC& rho, NUMERIC& epsilon, const NUMERIC P) {
        rho = P - this->P_f + this->eps_f;
        epsilon = 0._num;
}

NUMERIC CausalEoS::min_P() {
    return 0._num;
}
NUMERIC CausalEoS::min_rho() {
    return 0._num;
}
NUMERIC CausalEoS::min_e() {
    return 0._num;
}


/* EffectiveBosonicEoS */
NUMERIC EffectiveBosonicEoS::get_P_from_rho(const NUMERIC rho_in, const NUMERIC epsilon) {
	// (note: if for some reason, this function causes problems, then it is likely better to implement a root-finding algorithm to find Phi_eff))
	// It is possible to compute P(rho) by inverting: rho = 2* sqrt(mu^2 + lambda*phi^2)*phi^2
	// and then using the energy density e, to compute P(e) using the EoS:

	// use the cubic formula for: 0 = a*y^2 + y^3 - b (see: https://www.wolframalpha.com/input?i=find+root+of+a*x%5E2+%2B+x%5E3+-+b ),
	// where y= phi^2 and a& be are defined as:
	NUMERIC a = std::pow(this->mu,2) / this->lambda;
	NUMERIC b = rho_in*rho_in / (4._num*this->lambda);

	NUMERIC inner_root = std::sqrt(27._num*b*b - 4._num*a*a*a*b); // helper variable
	NUMERIC long_term = std::cbrt( 3._num*std::sqrt(3.)*inner_root -2._num*a*a*a + 27._num*b ) ; // another helper variable
	// calc the effective Phi:
	NUMERIC y = ( long_term/std::cbrt(2._num) + std::cbrt(2._num)*a*a / long_term - a ) / 3._num;	// cubic formula
	NUMERIC Phi_eff = std::sqrt(y);
	// compute energy density from the effective Phi-field:
	NUMERIC e_tot = 2._num* std::pow(mu ,2)* std::pow(Phi_eff,2) + 1.5_num*lambda*std::pow(Phi_eff,4);

	// return the P_tablesure P(e) = P(e(Phi)) = P(e(Phi(rho))) => P(rho)
	return ( get_P_from_e(e_tot));
}

NUMERIC EffectiveBosonicEoS::get_P_from_e(const NUMERIC etot_in) {
	// etot is the total energy density of the fluid
	return ( (4._num/9._num)*this->rho0*std::pow( std::sqrt(1._num+ (3._num/4._num)*(etot_in/this->rho0) )-1._num, 2) );	// p = 4/9 * rho0 * ( sqrt(1 + 3/4 * rho/rho0) -1 )^2
}

NUMERIC EffectiveBosonicEoS::get_e_from_P(const NUMERIC P_in) {
	// etot is the total energy density of the fluid
	return ( 3._num*P_in + 4._num* std::sqrt( P_in * this->rho0) );	// positive root taken fron rho= 3*P +/- 4* sqrt(P*rho0) );	// p = 4/9 * rho0 * ( sqrt(1 + 3/4 * rho/rho0) -1 )^2
}

NUMERIC EffectiveBosonicEoS::dP_de(const NUMERIC etot_in) {
	// etot is actually the total energy density of the fluid
	return ( 1._num/3._num - std::pow(1._num+ (3._num/4._num)*(etot_in/this->rho0) , -0.5_num) / 3.0_num );
}

NUMERIC EffectiveBosonicEoS::dP_drho(const NUMERIC rho_in, const NUMERIC epsilon) {
	// is not implemented yet (only really viable to do this numerically, though)
	return 0.0_num;
}

void EffectiveBosonicEoS::callEOS(NUMERIC& myrho, NUMERIC& epsilon, const NUMERIC P) {
	// update restmass density (rho) and specific energy (epsilon) according to effective bosonic EOS
    NUMERIC e_tot = 3._num*P + 4._num* std::sqrt( P * this->rho0);	// positive root taken fron etot= 3*P +/- 4* sqrt(P*rho0)

	// we can use the effective Phi-field to compute rho and epsilon.
	// compute phi from e_tot using: etot = 2 mu^2 phi^2 + 3/2 lambda phi^4 (valid for all radii):
	// solve the equation using the quadratic formula and only take the positive roots since phi > 0 everywhere
	NUMERIC p = 4._num*std::pow(this->mu, 2) / (3._num * this->lambda);
	NUMERIC q = -2._num*e_tot / (3._num*this->lambda);
	// compute the effective phi:
	NUMERIC Phi_eff = std::sqrt( std::sqrt( std::pow(0.5_num*p,2) - q ) - 0.5_num*p );	// is always > 0

	// now compute rho = 2*sqrt( mu^2 + lambda phi^2) * phi^2
	myrho = 2._num*std::sqrt( std::pow(this->mu,2) + this->lambda*std::pow(Phi_eff,2) ) * std::pow(Phi_eff,2);

	// now use: etot = rho*(1+epsilon)
    epsilon = e_tot/myrho - 1._num;

	return;
}

NUMERIC EffectiveBosonicEoS::min_P() {
    return 0._num;
}
NUMERIC EffectiveBosonicEoS::min_rho() {
    return 0._num;
}
NUMERIC EffectiveBosonicEoS::min_e() {
    return 0._num;
}

NUMERIC EffectiveBosonicEoS::get_mu() {
	return this->mu;
}
NUMERIC EffectiveBosonicEoS::get_lambda() {
	return this->lambda;
}

/******************
 *  EoStable *
 ******************/


/* This function expects a filename to a table.
 * The table has to contain the rest mass density rho, P_tablesure P, energy density e
 * The indices of these values is given by the std::map such that
 *       indices["rho"] corresponds to the column of rho [1/fm^3]
 *       indices["e"] corresponds to the column of energy density [MeV/fm^3]
 *       indices["P"] corresponds to the column of P_tablesure [MeV/fm^3]
 *
 * Columns are counted starting from 0
 * Lines containing # are ignored
 */
bool EoStable::load_from_file(const std::string filename, std::map<std::string, int> indices) {

    const NUMERIC MeV_fm3_to_codeunits =	2.886376934e-6_num; // unit conversion from Mev/fm^3 to code units M_s*c^2/ (G*M_s/c^2)^3
    const NUMERIC neutron_mass = 939.565379_num;	 // nuclear density in MeV

    std::ifstream infile;
    infile.open(filename);
    if(!infile.is_open())
        return false;

    rho_table.clear(); P_table.clear(); e_table.clear();
    int max_index = 0;
    max_index = std::max( std::max(indices.at("rho"), indices.at("P")), indices.at("e"));

    std::cout << "loading EoS table " << filename << std::endl; // TODO: Check error that arises when this is commented out xD
    std::string line;
    while (std::getline(infile, line)) {
        NUMERIC temp;
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
 * | restmass density rho [1/fm^3] |  (skipped)  | energy density[MeV/fm^3]  |  P_tablesure[MeV/fm^3] |
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
 * If we are below or above the P_tablesures in P_table, 0. is returned
 *  otherwise simple linear interpolation is used to obtain rho, epsilon */
void EoStable::callEOS(NUMERIC& rho, NUMERIC& epsilon, const NUMERIC P) {

	// use P to scan the EOS table (iterate from the top first because we start at high P_tablesures):
	unsigned int table_len = P_table.size();
	NUMERIC e_tot_tmp = 0._num;


	if (P < P_table[0]) { // we are outside the validity range of the table. Return zero
		rho = 0._num;
		epsilon = 0._num;
		return;
	}

	// search the table from the smaller values on first
	for (unsigned int i = 1; i<table_len; i++) {
		if (P_table[i] > P) {
			// the correct value is between the i-1th index and the ith index:
			// interpolate linearily between them:
			rho = rho_table[i-1] + (rho_table[i] - rho_table[i-1]) / (P_table[i] - P_table[i-1]) * (P - P_table[i-1]);
			e_tot_tmp = e_table[i-1] + (e_table[i] - e_table[i-1]) / (P_table[i] - P_table[i-1]) * (P - P_table[i-1]);
			epsilon = e_tot_tmp/rho - 1._num;	// re-arrange to get epsilon. e=rho*(1+epsilon)
			return;
		}
	}
	// return 0. outside the validity
    rho = 0._num;
    epsilon = 0._num;
}

/* This expects as input rho and returns P (epsilon is ignored)
 * according to the tabulated EoS
 * If we are below or above the densities in rho, 0. is returned
 *  otherwise simple linear interpolation is used to obtain P from rho */
NUMERIC EoStable::get_P_from_rho(const NUMERIC rho, const NUMERIC epsilon) {
    unsigned int table_len = rho_table.size();

    // if we are below the table return 0.
    if (rho < rho_table[0]) {
        return 0._num;
	}

	// search the table for the correct value
	for (unsigned int i = 1; i<table_len; i++) {
		if (rho_table[i] > rho) {
			// the correct value is between the ith index and the i+1th index:
			// interpolate linearily between them:
			return P_table[i-1] + (P_table[i] - P_table[i-1]) / (rho_table[i] - rho_table[i-1]) * (rho - rho_table[i-1]);
		}
	}
    // if no value was found extrapolate linearly
	return P_table[table_len-2] + (P_table[table_len-1] - P_table[table_len-2]) / (rho_table[table_len-1] - rho_table[table_len-2]) * (rho - rho_table[table_len-2]);
}

// obtain P from an e_tot input
NUMERIC EoStable::get_P_from_e(const NUMERIC e) {
	unsigned int table_len = e_table.size();

    // if we are below the table interpolate between (0.,0.) and (rho[0], P[0])
    if (e <= e_table[0]) {
        return 0._num + (P_table[0]-0._num) / (e_table[0]-0._num) * (e-0._num);
    }

	// search the table for the correct value
	for (unsigned int i = 1; i<table_len; i++) {
		if (e_table[i] > e) {
			// the correct value is between the ith index and the i+1th index:
			// interpolate linearily between them:
			return P_table[i-1] + (P_table[i] - P_table[i-1]) / (e_table[i] - e_table[i-1]) * (e - e_table[i-1]);
		}
	}
    // if no value was found extrapolate linearly
	return P_table[table_len-2] + (P_table[table_len-1] - P_table[table_len-2]) / (e_table[table_len-1] - e_table[table_len-2]) * (e- e_table[table_len-2]);
}

// obtain P from an e-tot input
NUMERIC EoStable::get_e_from_P(const NUMERIC P) {
	unsigned int table_len = P_table.size();

    // if we are below the table interpolate between (0.,0.) and (rho[0], P[0])
    if (P <= P_table[0]) {
        return 0._num + (e_table[0]-0._num) / (P_table[0]-0._num) * (P-0._num);
    }

	// search the table for the correct value
	for (unsigned int i = 1; i<table_len; i++) {
		if (P_table[i] > P) {
			// the correct value is between the ith index and the i+1th index:
			// interpolate linearily between them:
			return e_table[i-1] + (e_table[i] - e_table[i-1]) / (P_table[i] - P_table[i-1]) * (P - P_table[i-1]);
		}
	}
    // if no value was found extrapolate linearly
	return e_table[table_len-2] + (e_table[table_len-1] - e_table[table_len-2]) / (P_table[table_len-1] - P_table[table_len-2]) * (P - P_table[table_len-2]);
}

// obtain dP/drho from a rho input
NUMERIC EoStable::dP_drho(const NUMERIC rho, const NUMERIC epsilon) {
	// search the table for the correct value:
	unsigned int table_len = rho_table.size();

    //assert(rho_in < rho[table_len-2]); // out of range will return 0.
    if(rho < rho_table[1]) {
        NUMERIC dP1 = (P_table[2]-P_table[1])/(rho_table[2]-rho_table[1])/2._num + (P_table[1] - P_table[0])/(rho_table[1]-rho_table[0])/2._num;
        return 0._num + dP1 / (rho_table[1] - 0._num) * (rho - 0._num);
    }

	for (unsigned int i = 2; i<table_len; i++) {
		// scan for the first matching rho
		if (rho_table[i] > rho) {
            // estimate derivative at point i-1 and i
            NUMERIC dP1 = (P_table[i]-P_table[i-1])/(rho_table[i]-rho_table[i-1])/2._num + (P_table[i-1] - P_table[i-2])/(rho_table[i-1]-rho_table[i-2])/2._num;
            NUMERIC dP2 = (P_table[i+1]-P_table[i])/(rho_table[i+1]-rho_table[i])/2._num + (P_table[i] - P_table[i-1])/(rho_table[i]-rho_table[i-1])/2._num;
            return dP1 + (dP2-dP1)/(rho_table[i]- rho_table[i-1]) * (rho - rho_table[i-1]);
		}
	}
    // if no value was found return 0.
    return 0._num;
}

/* This expects as input rho, epsilon and returns dP/de
 * according to the tabulated EoS
 * If we are below or above the energies in e, 0. is returned
 *  otherwise simple linear interpolation is used to obtain P from e */
NUMERIC EoStable::dP_de(const NUMERIC e) {
    // search the table for the correct value:
    unsigned int table_len = e_table.size();

    if(e < e_table[1])
        return 0._num;

    for (unsigned int i = 2; i<table_len; i++) {
        // scan for the first matching e
        if (e_table[i] > e) {
            // estimate derivative at point i-1 and i
            NUMERIC dP1 = (P_table[i]-P_table[i-1])/(e_table[i]-e_table[i-1])/2._num + (P_table[i-1] - P_table[i-2])/(e_table[i-1]-e_table[i-2])/2._num;
            NUMERIC dP2 = (P_table[i+1]-P_table[i])/(e_table[i+1]-e_table[i])/2._num + (P_table[i] - P_table[i-1])/(e_table[i]-e_table[i-1])/2._num;
            return dP1 + (dP2-dP1)/(e_table[i]- e_table[i-1]) * (e - e_table[i-1]);
        }
    }
    return 0.;
}

NUMERIC EoStable::min_P() {
    return this->P_table.at(0);
}
NUMERIC EoStable::min_rho() {
    return this->rho_table.at(0);
}
NUMERIC EoStable::min_e() {
    return this->e_table.at(0);
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
