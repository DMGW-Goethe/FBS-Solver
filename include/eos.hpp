#ifndef EOS_HPP
#define EOS_HPP

#include <vector>
#include <cmath>
#include <fstream>	// filestream for file input
#include <sstream>	// string stream used for reading the EOS table
#include <iostream>
#include <stdexcept> // for std::runtime_error
#include <map>
#include "vector.hpp" // for NUMERIC

// a class to contain tabulated EOS
// reads in a EOS table in the constructor.
// Then it serves as a look-up table for the EOS (p=p(rho,epsilon,...)) to use in the integrator.

/* Equation of State
 *  an abstract class to model what an equation of state should contain
 *
 *  min_P, min_rho describe minimal P, minimal rho where the EoS is valid
 * */
class EquationOfState
{
public:
    virtual NUMERIC get_P_from_rho(const NUMERIC rho, const NUMERIC epsilon)  = 0;
	virtual NUMERIC get_P_from_e(const NUMERIC e)  = 0;
	virtual NUMERIC get_e_from_P(const NUMERIC P)  = 0;
    virtual NUMERIC dP_drho(const NUMERIC rho, NUMERIC epsilon) = 0;
	virtual NUMERIC dP_de(const NUMERIC e) = 0;

    /* This gives the derivative dP/de depending on the matter density rho and internal energy epsilon */
    virtual NUMERIC dP_de(const NUMERIC rho, NUMERIC epsilon)
	{  return dP_de(rho*(1._num + epsilon));  }

    /* This gives the matter density rho and internal energy epsilon depending on the pressure P */
	virtual void callEOS(NUMERIC& myrho, NUMERIC& epsilon, const NUMERIC P) = 0;

	/* Functions giving minimal P and minimal rho for the EoS to be valid */
    virtual NUMERIC min_P() = 0;		// minimal value of pressure
    virtual NUMERIC min_rho() = 0;	// minimal value of restmass density
	virtual NUMERIC min_e() = 0;	// minimal value of total energy density e:=rho(1+epsilon)

};

/* PolytropicEoS
 * a class modeling a Polytropic equation of state
 * with the parameters Kappa, Gamma
 * such that
 *  P = Kappa * rho^{Gamma}
 * */
class PolytropicEoS : public EquationOfState
{
protected:
    NUMERIC kappa, Gamma;

public:
    PolytropicEoS(const NUMERIC kappa=100._num, const NUMERIC Gamma =2._num) : kappa(kappa), Gamma(Gamma) {}

    /* This gives the pressure P depending on the matter density rho. The internal energy density is ignored */
    NUMERIC get_P_from_rho(const NUMERIC rho_in, const NUMERIC epsilon);
	NUMERIC get_P_from_e(const NUMERIC etot_in);
	NUMERIC get_e_from_P(const NUMERIC P_in);

	NUMERIC dP_drho(const NUMERIC rho, NUMERIC epsilon);

    /* This gives the derivative dP/de depending on the matter density rho. The internal energy density is ignored */
    NUMERIC dP_de(const NUMERIC e);
    /* This gives the matter density rho and internal energy density depending on the pressure P */
	void callEOS(NUMERIC& myrho, NUMERIC& epsilon, const NUMERIC P);

    /* Functions giving minimal P and minimal rho. For the polytrope both are 0. */
    NUMERIC min_P();
    NUMERIC min_rho();
	NUMERIC min_e();

};


/* CausalEoS
 * a class modeling a causal equation of state
 * with parameters eps_f, P_F,
 * such that
 *  P = P_f + rho(1+eps) - eps_f
 * */
class CausalEoS : public EquationOfState
{
protected:
    NUMERIC eps_f, P_f;

public:
    CausalEoS(const NUMERIC eps_f, const NUMERIC P_f =0._num) : eps_f(eps_f), P_f(P_f) {}

    /* This gives the pressure P depending on the matter density rho and internal energy epsilon */
    NUMERIC get_P_from_rho(const NUMERIC rho_in, const NUMERIC epsilon);
	NUMERIC get_P_from_e(const NUMERIC etot_in);
	NUMERIC get_e_from_P(const NUMERIC P_in);

    /* This gives the derivative dP/de depending on the matter density rho and internal energy epsilon */
    NUMERIC dP_de(const NUMERIC e);
	NUMERIC dP_drho(const NUMERIC rho, NUMERIC epsilon);

    /* This gives the matter density rho and internal energy epsilon depending on the pressure P */
	void callEOS(NUMERIC& myrho, NUMERIC& epsilon, const NUMERIC P);

    /* Functions giving minimal P and minimal rho. For the causal EoS both are 0. */
    NUMERIC min_P();
    NUMERIC min_rho();
	NUMERIC min_e();
};


/* a class modeling the effective equation of state for a bosonic condensate of self-interacting bosons */
class EffectiveBosonicEoS : public EquationOfState
{
protected:
    NUMERIC rho0;	// parameter computed from boson mass and self-interaction-parameter, corresponds to total energy density of the boson-fluid
					// rho0 = mu^4 / ( 2 * lambda ) in our normalization-convention for Phi and lambda (different convention to the two papers below:)
					// compare to: PHYSICAL REVIEW LETTERS VOLUME 57, Number 20 17 NOVEMBER 1986 Boson Stars: Gravitational Equilibria of Self-Gravitating scalar fields
					// and: Tidal deformability of dark matter admixed neutron stars PHYSICAL REVIEW D 105, 123010 (2022)
	NUMERIC mu, lambda;	// make the variables part of the class explicitly
public:
    EffectiveBosonicEoS(const NUMERIC mu=1._num, const NUMERIC lambda =1._num) : rho0(std::pow(mu,4) / (2._num* lambda)), mu(mu), lambda(lambda) {}

    NUMERIC get_P_from_rho(const NUMERIC rho_in, const NUMERIC epsilon);
	NUMERIC get_P_from_e(const NUMERIC etot_in);
	NUMERIC get_e_from_P(const NUMERIC P_in);
	void callEOS(NUMERIC& myrho, NUMERIC& epsilon, const NUMERIC P);
    NUMERIC dP_drho(const NUMERIC rho, NUMERIC epsilon);
	NUMERIC dP_de(const NUMERIC etot);

	NUMERIC get_mu();
	NUMERIC get_lambda();

    NUMERIC min_P();
    NUMERIC min_rho();
	NUMERIC min_e();
};


/* EoStable
 * a class representing a tabulated equation of state
 * On construction, the tabulated EoS has to be loaded into the
 *  vectors rho, Pres, e_tot,   the mass density, Pressure, and energy density respectively
 * On call, the values are interpolated linearly from the table
 * */
class EoStable : public EquationOfState
{
protected:
	std::vector<NUMERIC> rho_table, P_table, e_table;

public:
    EoStable(const std::vector<NUMERIC>& rho_table, const std::vector<NUMERIC>& P_table, const std::vector<NUMERIC>& e_table)
        : rho_table(rho_table), P_table(P_table), e_table(e_table) {}

    /* Constructor expects link to file */
    EoStable(const std::string filename) : rho_table({}), P_table({}), e_table({})
        { if(! load_from_file(filename))
                throw std::runtime_error("File '" + filename + "' could not be loaded") ;  }

    /* This function loads an EoS from file, the first one with rigid column indices, the second one with arbitrary indices*/
    bool load_from_file(const std::string filename);
    bool load_from_file(const std::string filename, std::map<std::string, int> indices);

    /* This gives the pressure P depending on the matter density rho and internal energy epsilon by linear interpolation of the table*/
	void callEOS(NUMERIC& myrho, NUMERIC& epsilon, const NUMERIC P);
    NUMERIC dP_drho(const NUMERIC rho, NUMERIC epsilon);
	NUMERIC dP_de(const NUMERIC e);
    /* This gives the derivative dP/de depending on the matter density rho and internal energy epsilon by linear interpolation */
    /* This gives the matter density rho and internal energy epsilon depending on the pressure P by linear interpolation*/
    NUMERIC get_P_from_rho(const NUMERIC rho, const NUMERIC epsilon);
	NUMERIC get_P_from_e(const NUMERIC e);
	NUMERIC get_e_from_P(const NUMERIC P);

    /* Functions giving minimal P and minimal rho. These depend on the lowest values in the tables */
    NUMERIC min_P();
    NUMERIC min_rho();
	NUMERIC min_e();

};

#endif
