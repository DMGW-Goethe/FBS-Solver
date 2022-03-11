# How to use the code:

## Compile and run:

I recommend using an IDE such as VScode to compile an run. I personally use the cmake extension in VScode.
Alternatively, just compile the code using g++ directly:
compile:
g++ *.cpp
(this will compile all cpp files in the directory)
to run the code:
./a.out
(execute the a.out file, or run a.exe if you are on Windows)

## The custom functions:

In general: the code will integrate the system of coupled ODEs defined in ODE_system().
To do that, we use a Runge-Kutta-45-Fehlberg integrator. It is accurate to 5th order in the truncation error (taylor series cut off).
It has the advantage that it is more accurate than the normal Runge-Kutta scheme and it is easy to implement a variable step size.

In Shooting_integrator_save_intermediate() and Shooting_integrator_nosave_intermediate() we integrate the ODEs using Fehlberg but we also perform the shooting method with regard to omega (scalar field frequency). To find omega we use a Newton-Raphson root-finding method for the function Phi(omega, r->inf) = 0 (we impose a vanishing scalar field at r=infinity).
The difference between ...save_intermediate and ...nosave_intermediate is that save_intermediate requires an additional 2D array as an input. It will then save the result at each step in the array. This way, we can obtain radial profiles of all evolved quantities with respect to r. nosave_intermediate will only output the fermionic radius and the total (ADM) mass. If we are only interested in creating an MR relation, we should use the nosave_intermediate function because it does not fill an array with values each time step. If we want to obtain a radal profile for ONE star, we should use save_intermediate.

save_data_ode() will write the radial profile data obtained from Shooting_integrator_save_intermediate() into a txt file.
save_data_MR() will write one MR curve into a txt file.

## the custom classes:

There are two custom classes:
- vec5d is just a 5-dim vector class which holds double precision numbers. dfault vector operations are defined (e.g. addition, multiplication with scalar, scalar product, ...)
- eostable is a class to hold the data of one tabulated EOS for restmass density, pressure and energy density. When reading the data from a table, the values are directly converted to code units. The table is loaded in when initializing the object. A member function is defined to obtain epsilon and rho from a given P (pressure) input.

## To use the functions as a black box:

- To integrate one star without radial profile output:

call: Shooting_integrator_nosave_intermediate(R_fermi, M_total, inits, myEOS, mu, lambda)
R_fermi and M_total will be changed by the function by pointer.
inits is a vec5d which holds the initial values at r=0 for  (a, alpha, Phi, Psi, P(rho_c)) in that order.
myEOS is an instance of the eostable class.
mu and lambda are mass (mu) and self-interaction parameter (lambda) of the dark matter.

- To integrate one star WITH radial profile output:
same as above but this time call: Shooting_integrator_save_intermediate(ODEdata, Nsteps, R_fermi, M_total, inits, myEOS, mu, lambda).
myarray is a 2D array of the form double ODEdata[N_steps+1][6]. We integrate a number of N_steps (array is 1 longer because we also save the init values) at each step i we then have inside of the ODEdata-array: ODEdata[i] = (r, a, alpha, Phi, Psi, P).

- To create an MR curve:

call: Shooting_integrator_nosave_intermediate(R_fermi, M_total, inits, myEOS, mu, lambda)
where: inits = (a, alpha, Phi, Psi, P(rho_c))
inside of a for-loop and change rho_c (central density) each iteration step. Then save the R and M values in a suitable array.
P(rho_c) can either be obtained using a polytropic EOS or using the eostable class member function get_P_from_rho(rho_c).