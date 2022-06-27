import numpy as np
from scipy.optimize import fsolve # used to calculate the radius of the NS matter (only the baryonic part)

import integrator
import eos_importer
import helper_functions
import eigenvalue_finder
import parameters

import importlib
importlib.reload(integrator)
importlib.reload(eos_importer)
importlib.reload(helper_functions)
importlib.reload(eigenvalue_finder)
importlib.reload(parameters)

class FieldConfiguration:
    def __init__(self, phic, rhoc, omega = -1) -> None:
        
        self.EoS = parameters.EoS_v

        self.phic = phic
        self.rhoc = rhoc
        
        self.omega = omega
        self.omegaMin = -1
        self.omegaMax = -1
        self.performedShooting = False

        self.totGravitationalMass = -1
        
        self.fNumber = -1
        self.bNumber = -1
        self.fRadius = -1
        
        self.BtoFratio = -1
        self.FtoBratio = -1

        self.convergenceRadius = -1
        self.convergenceId = -1

        self.terminator = "none"
        self.sol = "none"

    def radVals(self):
        return self.sol['t']
    
    def aVals(self):
        return self.sol['y'][0]
    
    def alphaVals(self):
        return self.sol['y'][1]

    def phiVals(self):
        return self.sol['y'][2]

    def psiVals(self):
        return self.sol['y'][3]

    def pVals(self):
        return self.sol['y'][4]
    
    def calculateProtoSolution(self, omega = -1):
        if omega == -1:
            self.sol = integrator.getProtoSolution(self.phic, self.rhoc, self.omega)
        
        else:
            self.sol = integrator.getProtoSolution(self.phic, self.rhoc, omega)

        return self.sol

    def findSolution(self, omegaGuess = 1, omegaMinGuess = 0, omegaMaxGuess = 2):
        if self.performedShooting:
            return self.sol
        
        if self.omegaMin != -1:
            omegaMinGuess = self.omegaMin
        if self.omegaMax != -1:
            omegaMaxGuess = self.omegaMax

        self.omega, self.omegaMin, self.omegaMax = eigenvalue_finder.binarySearch(self.phic, self.rhoc, omegaMinGuess, omegaMaxGuess, eigenvalue_finder.terminateMass)
        self.calculateProtoSolution()

    def getFermionRadius(self):
        if(self.rhoc == 0):
            return 0

        if self.fRadius != -1:
            return self.fRadius

        pressure = self.pVals()
        radii = self.radVals()
        minPressure = self.EoS.minPressure


        for i in range(len(pressure)):
            if(pressure[i] <= minPressure):
                return radii[i] * 1.477
        # minPressure = 10e-10

        #interpolated = lambda r : np.interp(r, radii, P) - minPressure#
        # TODO: the following line sometimes unstable. I am also not happy with the xtol and would like to have it smaller
        solution, infodict, ier, mesg = fsolve(lambda r : np.interp(r, radii * 1.477, pressure) - minPressure, 1, xtol = 1e-6, full_output = True)
        
        if ier != 1:    
            solution, infodict, ier, mesg = fsolve(lambda r : np.interp(r, radii * 1.477, pressure) - minPressure, 10, xtol = 1e-6, full_output = True)
            if ier != 1:
                print("Error: (getFermionRadius) can not calculate fRadius", mesg)
                return -1

        self.fRadius = solution[0]
        return self.fRadius

    def getGravitationalMass(self):
        if self.totGravitationalMass != -1:
            return self.totGravitationalMass

        rVals = self.radVals()
        aVals = self.aVals()
        amtVals = len(rVals)

        massVals = [rVals[i] / 2 * (1 - 1 / aVals[i]**2) for i in range(len(rVals))]
        derivVals = []

        for i in range(amtVals - 1):
            deriv = (massVals[i + 1] - massVals[i])/(rVals[i+1] - rVals[i])
            derivVals.append(deriv)

        minInd, minVal = helper_functions.localMinima(derivVals)

        self.totGravitationalMass = massVals[minInd]
        self.convergenceRadius = rVals[minInd]
        self.convergenceId = minInd

        return massVals[minInd]

    def getBaryonNumber(self, maxIndex = -1):
        if self.bNumber != -1:
            return self.bNumber

        radii = self.radVals()[0:maxIndex]
        a = self.aVals()[0:maxIndex]
        alpha = self.alphaVals()[0:maxIndex]
        phi = self.phiVals()[0:maxIndex]

        arr1 = [4 * np.pi * self.omega * a[i] * phi[i]**2 * radii[i]**2 / alpha[i] for i in range(len(radii))]

        numBaryons = np.trapz(arr1, x = radii)

        self.bNumber = numBaryons
        return numBaryons
    
    def getFermionNumber(self):
        if self.rhoc == 0:
            return 0
        if self.fNumber != -1:
            return self.fNumber

        radii = self.radVals()
        a = self.aVals()
        pressure = self.pVals()
        
        radius = self.getFermionRadius()
        if radius == -1:
            print("Error: can not find fermion radius")
            return -1

        maxFermionIndex = 0
        
        for i in range(len(radii)):
            if radii[i] * 1.477 > radius:
                maxFermionIndex = i
                break

        restmassrho = [self.EoS.get_restMassDensity(p) for p in pressure]
        # restmassrho = [EoS_polytrope_rho(p) for p in pressure]

        arr2 = [4 * np.pi * a[i] * restmassrho[i] * radii[i]**2 for i in range(len(radii))]

        numFermions = np.trapz(arr2[0:maxFermionIndex + 1], x = radii[0:maxFermionIndex + 1])
        return numFermions

# parameters.EoS_v = eos_importer.EoS("DD2")
# test2 = FieldConfiguration(0.04, 0.002)

# test2.findSolution()
# print(test2.getGravitationalMass())


# test.findSolution()

# import matplotlib.pyplot as plt

# plt.plot(test.radVals(), test.pVals())
# plt.yscale('log')
# plt.show()

# print(test.pVals()[-1])
# print(test.getFermionRadius())





# # this is the >total< gravitational mass (fermions + bosons)
# def getGravitationalMass(solution):
#     aend = solution['y'][0][-1]
#     alphaend = solution['y'][1][-1]
#     rend = solution['t'][-1]

#     return rend / 2 * (1 - 1 / aend**2)
