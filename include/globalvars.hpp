#ifndef GLOBALVARS_HPP
#define GLOBALVARS_HPP

// all global variables are contained here.
// -----------------------------------------------------
// by changing this, one can chose the wanted floating-point precision.
// "out-of the box" supported data types are 'double' and 'long double'
#define NUMERIC double
// user-defined literal to automatically convert every floating-point constant to NUMERIC type
inline NUMERIC operator "" _num(long double d) {return (NUMERIC)d;}

// parameters used in the "fbs"-class. Used in the "fps.cpp"-file
#define PHI_converged 1e-4_num // compared to phi / phi_0, when the bosonic component has converged sufficiently
#define INT_converged 1e-7_num // integration converged
#define M_T_converged 1e-15_num // leftover from previous attempts to characterize convergence

// parameters related to the "nsmodel"-class. used in the "nsmodel.cpp"-file and to functions related to the NS part of a FBS
#define R_INIT 1e-10_num    // the initial integration radius
#define R_MAX 500._num      // the general maximum integration radius (might be increased if not sufficient)
#define P_ns_min 1e-15_num  // the minimum pressure for the "boundary" of the NS

#endif