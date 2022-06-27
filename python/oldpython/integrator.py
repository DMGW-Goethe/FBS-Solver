from matplotlib.pyplot import getp
import numpy as np
from scipy.integrate import solve_ivp # I use solve_ivp instead of odeint. odeint might be faster, but solve_ivp allows to specify the used integration method
from scipy.optimize import fsolve
from scipy.integrate import odeint

import eos_importer
import parameters

import importlib
importlib.reload(eos_importer)
importlib.reload(parameters)


# mu = 1 # the mass of the bosonic field
# lam = 0 # the self-interaction parameter in the quartic coupling

# this functions returns df/dr (we assume spherical symmetry in all components)
def fun(r, y):
    # using the same order and notation as in the paper
    a = y[0] # this the radial factor in the metric
    alpha = y[1] # this the time factor in the metric
    phi = y[2] # this the field value of the bosonic field 
    psi = y[3] # this just the radial derivative of phi (= dphi/dr)
    pressure = y[4] # this the pressure
    omega = y[5] # I am cheating a bit by passing omega as one of the array components (not sure how to do it otherwise)

    rho = 0

    if(pressure < 0): # make sure that the pressure vanishes after we cross the surface (which is defined at P = 0). TODO: I actually define the surface with P < minpressure. Maybe I should change this criterion to if(p < minpressure) ?
        pressure = 0
    else:
        rho = parameters.EoS_v.get_restMassDensity(pressure)
        # rho = EoS_v.get_totalEnergyDensity(pressure)
    epsilon = 0
    if(rho != 0):
        epsilon = parameters.EoS_v.get_totalInternalEnergy(pressure)
        # epsilon = 0

    dadr = a / 2 * ((1 - a*a) / r + 4 * np.pi * r * ((omega**2/alpha**2 + parameters.mu**2 + parameters.lam/2 * phi**2) * a**2 * phi**2 + psi**2 + 2 * a**2 * rho * (1 + epsilon)))
    dalphadr = alpha / 2 * ((a*a - 1) / r + 4 * np.pi * r * ((omega**2/alpha**2 - parameters.mu**2 - parameters.lam/2 * phi**2) * a**2 * phi**2 + psi**2 + 2 * a**2 * pressure))
    dphidr = psi
    dpsidr = -(1 + a**2 - 4 * np.pi * r**2 * a**2 * (parameters.mu**2 * phi**2 + parameters.lam/2 * phi**4 + rho * (1 + epsilon) - pressure)) * psi / r - (omega**2 / alpha**2 - parameters.mu**2 - parameters.lam * phi**2) * a**2 * phi
    dPdr = -(rho * (1 + epsilon) + pressure) * dalphadr / alpha

    return [dadr, dalphadr, dphidr, dpsidr, dPdr, 0]

# this returns a solution to the system of ODEs in terms of phic, rhoc and omega. This is not yet a physical solution though. We have to choose omega so that phic vanishes at infinity
def getProtoSolution(phic, rhoc, omega):
    # initial conditions
    # a = 0, alpha = 0, phi = phi_c, psi = 0, p = K * rho_c^gamma, rho = rho_c (all quantities evaluated at r = 0)
    ac = 1
    alphac = 1
    psic = 0

    Pc = fsolve(lambda x : parameters.EoS_v.get_restMassDensity(x) - rhoc, [0.0])[0]
    # print(Pc)
    # Pc = EoS_polytrope_pressure(rhoc)

    y0 = [ac, alphac, phic, psic, Pc, omega]
    # rvals = np.linspace(r0, rmax, 100000)

    # ---> I tried to stop integrating when psi explodes too much, but this seems to actually slow down everything because solve_ivp then does a root finding algorithm to find out where exactly psi = xxx :'(
    # def event(t, y):
    #     return np.abs(y[3] - 100)

    # event.terminal = True

    return solve_ivp(fun, [parameters.minIntegrationRadius, parameters.maxIntegrationRadius], y0, atol = 1e-23, rtol = 1e-17, events = None)
    # return odeint(func = fun, y0 = y0, t = np.linspace(parameters.minIntegrationRadius, parameters.maxIntegrationRadius, 1000), atol = 1e-20, rtol = 1e-10, full_output = True, tfirst = True)
    

if __name__ == "__main__":

    test = getProtoSolution(0.04, 0.002, 1.2)
    # print(test)

    print("finihsed")