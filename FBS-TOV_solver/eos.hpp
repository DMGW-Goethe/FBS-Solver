#ifndef EOS_HPP
#define EOS_HPP

#include <vector>
#include <cmath>
#include <fstream>	// filestream for file input
#include <sstream>	// string stream used for reading the EOS table
#include <iostream>

// a class to contain tabulated EOS
// reads in a EOS table in the constructor.
// Then it serves as a look-up table for the EOS (p=p(rho,epsilon,...)) to use in the integrator.

class EquationOfState
{
public:
    virtual double get_P_from_rho(const double rho_in)  = 0;

	virtual void callEOS(double& myrho, double& epsilon, const double P) = 0;

};

class PolytropicEoS : public EquationOfState
{
protected:
    double kappa, Gamma;

public:
    PolytropicEoS(const double kappa=100., const double Gamma =2.);

    double get_P_from_rho(const double rho_in);
	void callEOS(double& myrho, double& epsilon, const double P);

};

class EoStable : public EquationOfState
{
public:
	// constructors:
    EoStable(const std::string filename);   // chose which EOS to load in using the filename and fill the std::vectors

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
