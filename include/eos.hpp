#ifndef EOS_HPP
#define EOS_HPP

#include <vector>
#include <cmath>
#include <fstream>	// filestream for file input
#include <sstream>	// string stream used for reading the EOS table
#include <iostream>
#include <stdexcept> // for std::runtime_error

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
    /* This gives the pressure P depending on the matter density rho and internal energy epsilon */
    virtual double get_P_from_rho(const double rho, const double epsilon)  = 0;

    /* This gives the derivative dP/drho depending on the matter density rho and internal energy epsilon */
    virtual double dP_drho(const double rho, double epsilon) = 0;

    /* This gives the matter density rho and internal energy epsilon depending on the pressure P */
	virtual void callEOS(double& myrho, double& epsilon, const double P) = 0;

    /* Functions giving minimal P and minimal rho for the EoS to be valid */
    virtual double min_P() = 0;
    virtual double min_rho() = 0;

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
    double kappa, Gamma;

public:
    PolytropicEoS(const double kappa=100., const double Gamma =2.) : kappa(kappa), Gamma(Gamma) {}

    /* This gives the pressure P depending on the matter density rho. The internal energy density is ignored */
    double get_P_from_rho(const double rho_in, const double epsilon);
    /* This gives the derivative dP/drho depending on the matter density rho. The internal energy density is ignored */
    double dP_drho(const double rho, double epsilon);
    /* This gives the matter density rho and internal energy density depending on the pressure P */
	void callEOS(double& myrho, double& epsilon, const double P);

    /* Functions giving minimal P and minimal rho. For the polytrope both are 0. */
    double min_P();
    double min_rho();

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
    double eps_f, P_f;

public:
    CausalEoS(const double eps_f, const double P_f =0.) : eps_f(eps_f), P_f(P_f) {}

    /* This gives the pressure P depending on the matter density rho and internal energy epsilon */
    double get_P_from_rho(const double rho_in, const double epsilon);

    /* This gives the derivative dP/drho depending on the matter density rho and internal energy epsilon */
    double dP_drho(const double rho, double epsilon);

    /* This gives the matter density rho and internal energy epsilon depending on the pressure P */
	void callEOS(double& myrho, double& epsilon, const double P);

    /* Functions giving minimal P and minimal rho. For the causal EoS both are 0. */
    double min_P();
    double min_rho();
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
	std::vector<double> rho_table, P_table, e_table;

public:
    /* Constructor expects link to file TODO: add constructor with lists */
    EoStable(const std::string filename);

    /* This gives the pressure P depending on the matter density rho and internal energy epsilon by linear interpolation of the table*/
	void callEOS(double& myrho, double& epsilon, const double P);
    /* This gives the derivative dP/drho depending on the matter density rho and internal energy epsilon by linear interpolation */
    double dP_drho(const double rho, double epsilon);
    /* This gives the matter density rho and internal energy epsilon depending on the pressure P by linear interpolation*/
    double get_P_from_rho(const double rho_in, const double epsilon);

    /* Functions giving minimal P and minimal rho. These depend on the lowest values in the tables */
    double min_P();
    double min_rho();

};

#endif
