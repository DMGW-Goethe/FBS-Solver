#ifndef EOS_HPP
#define EOS_HPP

#include <vector>
#include <cmath>
#include <fstream>	// filestream for file input
#include <sstream>	// string stream used for reading the EOS table
#include <iostream>
#include <cassert>	// for the assert() command

// a class to contain tabulated EOS
// reads in a EOS table in the constructor.
// Then it serves as a look-up table for the EOS (p=p(rho,epsilon,...)) to use in the integrator.

/* an abstract class to model what an equation of state should contain
 * */
class EquationOfState
{
public:
    virtual double get_P_from_rho(const double rho_in)  = 0;

	virtual void callEOS(double& myrho, double& epsilon, const double P) = 0;

};

/* a class modeling a Polytropic equation of state */
class PolytropicEoS : public EquationOfState
{
protected:
    double kappa, Gamma;

public:
    PolytropicEoS(const double kappa=100., const double Gamma =2.) : kappa(kappa), Gamma(Gamma) {}

    double get_P_from_rho(const double rho_in);
	void callEOS(double& myrho, double& epsilon, const double P);

};

/* a class modeling a tabulated equation of state */
class EoStable : public EquationOfState
{
public:
    EoStable(const std::string filename);
	// call EOS lookup-table (including interpolation):
	void callEOS(double& myrho, double& epsilon, const double P);
    double get_P_from_rho(const double rho_in);

private:
	// vectors to hold the tabulated EOS values:
	std::vector<double> rho;	// restmass density
	std::vector<double> Pres;	// pressure
	std::vector<double> e_tot;	// total energy density e:=rho(1+epsilon)
	// helper variables:
	unsigned current_index;		// saves position in look up table
	std::string eos_table_name; // name of EOS/table to use (filename of txt file)
};

#endif
