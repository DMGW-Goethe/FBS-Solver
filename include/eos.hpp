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

/* an abstract class to model what an equation of state should contain
 * */
class EquationOfState
{
public:
    virtual double get_P_from_rho(const double rho_in, const double epsilon)  = 0;
	virtual double get_P_from_etot(const double etot_in)  = 0;
	virtual double get_etot_from_P(const double P_in)  = 0;
    virtual double dP_drho(const double rho, double epsilon) = 0;
	virtual double dP_detot(const double etot) = 0;

	virtual void callEOS(double& myrho, double& epsilon, const double P) = 0;

    virtual double min_P() = 0;		// minimal value of pressure
    virtual double min_rho() = 0;	// minimal value of restmass density
	virtual double min_etot() = 0;	// minimal value of total energy density e:=rho(1+epsilon)

};

/* a class modeling a Polytropic equation of state */
class PolytropicEoS : public EquationOfState
{
protected:
    double kappa, Gamma;

public:
    PolytropicEoS(const double kappa=100., const double Gamma =2.) : kappa(kappa), Gamma(Gamma) {}

    double get_P_from_rho(const double rho_in, const double epsilon);
	double get_P_from_etot(const double etot_in);
	double get_etot_from_P(const double P_in);
	void callEOS(double& myrho, double& epsilon, const double P);
    double dP_drho(const double rho, double epsilon);
	double dP_detot(const double etot);

    double min_P();
    double min_rho();
	double min_etot();

};


/* a class modeling a Causal equation of state */
class CausalEoS : public EquationOfState
{
protected:
    double eps_f, P_f;

public:
    CausalEoS(const double eps_f, const double P_f =0.) : eps_f(eps_f), P_f(P_f) {}

    double get_P_from_rho(const double rho_in, const double epsilon);
    double dP_drho(const double rho, double epsilon);
	void callEOS(double& myrho, double& epsilon, const double P);

    double min_P();
    double min_rho();
	double min_etot();
};


/* a class modeling the effective equation of state for a bosonic condensate of self-interacting bosons */
class EffectiveBosonicEoS : public EquationOfState
{
protected:
    double rho0;	// parameter computed from boson mass and self-interaction-parameter, corresponds to total energy density of the boson-fluid
					// rho0 = mu^4 / ( 4 * lambda )

public:
    EffectiveBosonicEoS(const double mu=1.0, const double lambda =1.0) : rho0(std::pow(mu,4) / (4.* lambda)) {}

    double get_P_from_rho(const double rho_in, const double epsilon);
	double get_P_from_etot(const double etot_in);
	double get_etot_from_P(const double P_in);
	void callEOS(double& myrho, double& epsilon, const double P);
    double dP_drho(const double rho, double epsilon);
	double dP_detot(const double etot);

    double min_P();
    double min_rho();
	double min_etot();

};


/* a class modeling a tabulated equation of state */
class EoStable : public EquationOfState
{
public:
    EoStable(const std::string filename);
	// call EOS lookup-table (including interpolation):
	void callEOS(double& myrho, double& epsilon, const double P);
    double dP_drho(const double rho, double epsilon);
	double dP_detot(const double etot);
    double get_P_from_rho(const double rho_in, const double epsilon);
	double get_P_from_etot(const double etot_in);
	double get_etot_from_P(const double P_in);

    double min_P();
    double min_rho();
	double min_etot();

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
