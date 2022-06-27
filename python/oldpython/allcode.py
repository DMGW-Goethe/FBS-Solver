from typing import final
import numpy as np
from scipy.integrate import solve_ivp # I use solve_ivp instead of odeint. odeint might be faster, but solve_ivp allows to specify the used integration method
from scipy.optimize import fsolve # used to find the eigenvalue omega (i.e. carrying out the shooting method)
import matplotlib.pyplot as plt
import multiprocessing
from itertools import repeat

import sys

class EoS:
    neutron_mass = 939.565379 # MeV

    def __init__(self, name):
        self.name = name
        self.load()
    
    def load(self):
        self.raw_data = np.loadtxt(self.name + "_" + "eos.table").T

        # the following is in compose units
        self.baryonNumberDensity = self.raw_data[1] # fm^-3
        self.pressure = self.raw_data[3] # MeV fm^-3
        self.totalEnergyPerBaryon = self.raw_data[4] # MeV
        self.totalEnergyDensity = self.raw_data[4] # MeV fm^-3
        self.restMassDensity = self.raw_data[1] * self.neutron_mass # MeV fm^-3

        # add the 0 value for better convergence
        # TODO: add check if the data is already provided like this.
        # NO: this actually doesn't seem to work, because this will generate something that is not well behaved towards the origin
        # self.baryonNumberDensity = np.insert(self.baryonNumberDensity, 0, 0.0) # fm^-3
        # self.pressure = np.insert(self.pressure, 0, 0.0) # MeV fm^-3
        # self.totalEnergyPerBaryon = np.insert(self.totalEnergyPerBaryon, 0, 0.0) # MeV
        # self.totalEnergyDensity = np.insert(self.totalEnergyDensity, 0, 0.0) # MeV fm^-3
        # self.restMassDensity = np.insert(self.restMassDensity, 0, 0.0) # MeV fm^-3

        # the following is in code units (c = G = Msun = 1). Code units are depicted with cu
        self.pressure_cu = self.pressure * 2.88e-6
        # self.totalEnergyPerBaryon_cu = self.raw_data[4]
        self.totalEnergyDensity_cu = self.totalEnergyDensity * 2.88e-6
        self.restMassDensity_cu = self.restMassDensity * 2.886e-6
        self.amtEntries = len(self.pressure)

        self.totalEnergyDensity_vs_pressure_cu = lambda x : np.interp(x, self.pressure_cu, self.totalEnergyDensity_cu)
        self.restMassDensity_vs_pressure_cu = lambda x : np.interp(x, self.pressure_cu, self.restMassDensity_cu)

        self.minPressure = min(self.pressure_cu) # TODO: think of something better on where to take the fermion radius. Maybe I could interpolate starting from rho = 0 -> p = 0?
        # self.minPressure = 0

    def get_restMassDensity(self, pressure):        
        return self.restMassDensity_vs_pressure_cu(pressure)
        
    def get_totalEnergyDensity(self, pressure):        
        return self.totalEnergyDensity_vs_pressure_cu(pressure)
    
    def get_totalInternalEnergy(self, pressure):
        totEnergy = self.get_totalEnergyDensity(pressure)
        restMassDensity = self.get_restMassDensity(pressure)
        return totEnergy / restMassDensity - 1.0

# equation of state
# def EoS_polytrope_pressure(rho):
#     K = 100
#     gamma = 2
#     return K * rho**gamma

# def EoS_polytrope_rho(P):
#     K = 100
#     gamma = 2

#     if(P < 0.0):
#         return 0.0
#     return (P / K)**(1 / gamma)

class EoS_polytrope:
    def __init__(self):
        self.K = 100
        self.gamma = 2
        self.minPressure = 0
        pass

    def get_pressure(self, restMassDensity):
        return self.K * restMassDensity**self.gamma


    def get_restMassDensity(self, pressure):

        if(pressure < 0.0):
            return 0.0
        return (pressure / self.K)**(1 / self.gamma)
    
    def get_totalInternalEnergy(self, pressure):
        rho = self.get_restMassDensity(pressure)
        return self.K * rho**(self.gamma - 1) / (self.gamma - 1)

    def get_totalEnergyDensity(self, pressure):
        print("NOT IMPLEMENTED LOLOLLOL")
        return 0
    
import numpy as np
from scipy.integrate import solve_ivp # I use solve_ivp instead of odeint. odeint might be faster, but solve_ivp allows to specify the used integration method
from scipy.optimize import fsolve # used to find the eigenvalue omega (i.e. carrying out the shooting method)
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt

pi = np.pi
mu = 1 # the mass of the bosonic field
lam = 0 # the self-interaction parameter in the quartic coupling

# the min and max radius between which the ODEs are solved
minRadius = 1e-10
maxRadius = 200

# Here we load the DD2 equation of state. (The table with the data should be in the same directory as this notebook)
# EoS_v = EoS("DD2")
EoS_v = EoS_polytrope()

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
        rho = EoS_v.get_restMassDensity(pressure)
        # rho = EoS_v.get_totalEnergyDensity(pressure)
    epsilon = 0
    if(rho != 0):
        epsilon = EoS_v.get_totalInternalEnergy(pressure)
        # epsilon = 0

    dadr = a / 2 * ((1 - a*a) / r + 4 * pi * r * ((omega**2/alpha**2 + mu**2 + lam/2 * phi**2) * a**2 * phi**2 + psi**2 + 2 * a**2 * rho * (1 + epsilon)))
    dalphadr = alpha / 2 * ((a*a - 1) / r + 4 * pi * r * ((omega**2/alpha**2 - mu**2 - lam/2 * phi**2) * a**2 * phi**2 + psi**2 + 2 * a**2 * pressure))
    dphidr = psi
    dpsidr = -(1 + a**2 - 4 * pi * r**2 * a**2 * (mu**2 * phi**2 + lam/2 * phi**4 + rho * (1 + epsilon) - pressure)) * psi / r - (omega**2 / alpha**2 - mu**2 - lam * phi**2) * a**2 * phi
    dPdr = -(rho * (1 + epsilon) + pressure) * dalphadr / alpha

    return [dadr, dalphadr, dphidr, dpsidr, dPdr, 0]

# this returns a solution to the system of ODEs in terms of phic, rhoc and omega. This is not yet a physical solution though. We have to choose omega so that phic vanishes at infinity
def getProtoSolution(phic, rhoc, omega):
    # initial conditions
    # a = 0, alpha = 0, phi = phi_c, psi = 0, p = K * rho_c^gamma, rho = rho_c (all quantities evaluated at r = 0)
    ac = 1
    alphac = 1
    psic = 0

    Pc = fsolve(lambda x : EoS_v.get_restMassDensity(x) - rhoc, [0.0])[0]
    # print(Pc)
    # Pc = EoS_polytrope_pressure(rhoc)
    
    r0 = minRadius
    rmax = maxRadius

    y0 = [ac, alphac, phic, psic, Pc, omega]
    # rvals = np.linspace(r0, rmax, 100000)

    # ---> I tried to stop integrating when psi explodes too much, but this seems to actually slow down everything because solve_ivp then does a root finding algorithm to find out where exactly psi = xxx :'(
    # def event(t, y):
    #     return np.abs(y[3] - 100)

    # event.terminal = True

    return solve_ivp(fun, [r0, rmax], y0, atol = 1e-20, rtol = 1e-10, events = None)

# here we take a proto solution and turn it into a physical one by finding the value of omega such that phi goes to 0 for large radii
def getSolution(phic, rhoc):
    res = [2]
    
    # I first solve for omega to maximize the possible radius (the ode solver just stops at some radii if a problem occurs and this way we (somewhat) guarantee that we actually solved up to rmax)
    def root_rend(omega):
        guess = getProtoSolution(phic, rhoc, omega[0])
        rend = guess['t'][-1]
        phiend = guess['y'][2][-1]
        return abs(rend - maxRadius)

    res = fsolve(root_rend, [abs(res[0])])
    print(res)

    def get_aend(omega):
        sol = getProtoSolution(phic, rhoc, omega)
        aend = sol['y'][0][-1]
        return aend
    
    res = minimize_scalar(get_aend, bounds = [res[0] - 1, res[0] + 1], method = "bounded")

    print(res)
    return res

# returns value for the fermion (NS) radius
def getFermionRadius(solution):
    pressure = solution['y'][4]
    radii = solution['t']
    minPressure = EoS_v.minPressure
    # minPressure = 10e-10

    #interpolated = lambda r : np.interp(r, radii, P) - minPressure
    solution = fsolve(lambda r : np.interp(r, radii * 1.477, pressure) - minPressure, 10, xtol = 1e-22)
    return solution[0]

# this is the >total< gravitational mass (fermions + bosons)
def getGravitationalMass(solution):
    aend = solution['y'][0][-1]
    alphaend = solution['y'][1][-1]
    rend = solution['t'][-1]

    return rend / 2 * (1 - 1 / aend**2)

def getBaryonNumber(sol, omega, maxIndex = -1):
    radii = sol['t'][0:maxIndex]
    a = sol['y'][0][0:maxIndex]
    alpha = sol['y'][1][0:maxIndex]
    phi = sol['y'][2][0:maxIndex]

    arr1 = [4 * pi * omega * a[i] * phi[i]**2 * radii[i]**2 / alpha[i] for i in range(len(radii))]

    numBaryons = np.trapz(arr1, x = radii)
    return numBaryons

def getFermionNumber(sol):
    radii = sol['t']
    a = sol['y'][0]
    pressure = sol['y'][4]
    
    radius = getFermionRadius(sol)

    maxFermionIndex = 0
    
    for i in range(len(radii)):
        if radii[i] * 1.477 > radius:
            maxFermionIndex = i
            break

    restmassrho = [EoS_v.get_restMassDensity(p) for p in pressure]
    # restmassrho = [EoS_polytrope_rho(p) for p in pressure]

    arr2 = [4 * pi * a[i] * restmassrho[i] * radii[i]**2 for i in range(len(radii))]

    numFermions = np.trapz(arr2[0:maxFermionIndex + 1], x = radii[0:maxFermionIndex + 1])
    return numFermions


# function to find all roots in an array by iterating over all values and looking for changes in the sign of the value
def amtOfRoots(funArray):
    amtData = len(funArray)
    amtRoots = 0

    for i in range(amtData - 1):
        if np.sign(funArray[i]) != np.sign(funArray[i + 1]):
            amtRoots += 1
    
    return amtRoots

#returns the first local minima in the array given by vals
def localMinima(vals):
    for i in range(1, len(vals) - 1):
        if( vals[i - 1] > vals[i] and vals[i] < vals[i + 1] and vals[i] < 1e-5):
            return i, vals[i]
    return -1, -1

# binary search to find the omega value corresponding to a scalar field solution
# parameters: phic -> central scalar field value, omegaMin -> initial guess for minimal omega, omegaMax -> initial guess for the maximum omega
# the correct solution should be sandwiched by omegaMin and omegaMax in order for the algorithm to converge
#
# the function is structured in two parts: the first part will make sure that truly only the ground state (i.e. the state without any crossings with the r axis) is sandwiched by omegaMin and omegaMax
# this is done by decreasing omegaMax until only one root of the function phi(r) can be found within the integration interval
# the second part then does a binary search on the final value of phi: it checks wheter phi(r = rend) is positive or negative. Since the sign flip occurs when crossing the omega eigenvalue, we can use it to narrow down omega
def binarySearch(phic, rhoc, omegaMin, omegaMax, terminator = lambda x : False):
    if(phic == 0.0):
        return 0.0, 0.0, 0.0

    maxIteration = 100

    minSol = getProtoSolution(phic, rhoc, omegaMin)
    maxSol = getProtoSolution(phic, rhoc, omegaMax)

    minSolPhi = minSol['y'][2]
    maxSolPhi = maxSol['y'][2]

    # minSol should not cross phi = 0 at any point
    # maxSol should only cross once, so as a initial check we can decrease maxSol until only one intersection with the r axis is left
    if(amtOfRoots(minSolPhi) != 0):
        print("Error: the lower bound for omega should result in a distribution of phi with no roots at all")
        return 0.0, 0.0, 0.0

    # if the upper value of omega does not result in a root that is being sandwhiched, then we increase it until there is one
    while(amtOfRoots(maxSolPhi) == 0):
        omegaMin = omegaMax
        omegaMax += 0.1

        minSol = getProtoSolution(phic, rhoc, omegaMin)
        maxSol = getProtoSolution(phic, rhoc, omegaMax)

        minSolPhi = minSol['y'][2]
        maxSolPhi = maxSol['y'][2]

    # if(amtOfRoots(maxSolPhi) == 0):
    #     print("Error: the upper bound of omega should result in at least one root")
    #     return 0
    
    
    # Part1: decrease omegaMax until only the ground state is sandwiched by omegaMin and omegaMax
    # for this I define two temporary variables that will keep track on the upper and lower limit for omegaMax
    omegaMaxL = omegaMin
    omegaMaxU = omegaMax

    currAmtRoots = amtOfRoots(maxSolPhi) # this is how may roots are currently present for omega = omegaMax

    while True:
        # if we currently already have only one root, then we do not have to adjust omegaMax and can just skip this part
        if(currAmtRoots == 1):
            break
        
        # now we will test the value sad is in the middle of omegaMaxL and omegaMaxU and see how many roots are enclosed
        omegaMaxTemp = (omegaMaxL + omegaMaxU) / 2.0
        maxSolTemp = getProtoSolution(phic, rhoc, omegaMaxTemp)
        maxSolPhiTemp = maxSolTemp['y'][2]
        
        if(amtOfRoots(maxSolPhiTemp) == 0): # if we have zero roots, then the guess for omegaMax was too low and we have to increase it again. For this we can just set omegaMaxL as this value and continue
            omegaMaxL = omegaMaxTemp
            continue
        elif(amtOfRoots(maxSolPhiTemp) == 1): # if we have exactly one root then we are done
            omegaMaxU = omegaMaxTemp
            # print("finished")
            break
        else:
            omegaMaxU = omegaMaxTemp # this is triggered if we still have more than one root. in this case we set the new guess of omegaMaxU to the generated guess and continue
            continue
    

    # when we are done we can use the obtained guesses to also set omegaMin and omegaMax
    omegaMin = omegaMaxL
    omegaMax = omegaMaxU

    # print(omegaMin, omegaMax)

    # Part2: now we do a binary search on the last value of phi. (Note that we could also just keep the upper code running to do the binary search and also arrive at a value for omega. However, doing the binary search on the last value of phi will be faster because we don't have to iterate over the entire array to find all roots in every iteration)
    # initalize the values
    omegaMid = (omegaMin + omegaMax) / 2.0
    prevOmegaMid = omegaMid

    minSol = getProtoSolution(phic, rhoc, omegaMin)
    maxSol = getProtoSolution(phic, rhoc, omegaMax)

    for i in range(maxIteration):
        # the new guess will be generated from the value in the middle of omegaMin and omegaMax
        omegaMid = (omegaMin + omegaMax) / 2.0
        midSol = getProtoSolution(phic, rhoc, omegaMid)
        
        # save all signs for convenience
        minSign = np.sign(minSol['y'][2][-1])
        maxSign = np.sign(maxSol['y'][2][-1])
        midSign = np.sign(midSol['y'][2][-1])

        # if(i % 10 == 0):
        #     print("finished", i, "/", maxIteration)

        # check if we reached the desired accuracy
        # if(i != 0 and np.abs(prevOmegaMid - omegaMid) / omegaMid < 1e-30):
        #     # print(i)
        #     break
        if(terminator(midSol)):
            break

        prevOmegaMid = omegaMid

        # now we test which bound we have to adjust
        if(midSign == minSign):
            omegaMin = omegaMid
            minSol = midSol
            continue
        elif(midSign == maxSign):
            omegaMax = omegaMid
            maxSol = midSol
            continue
    
    return omegaMid, omegaMin, omegaMax

def terminateMass(sol):
    rVals = sol['t']
    aVals = sol['y'][0]
    amtVals = len(rVals)

    massVals = [rVals[i] / 2 * (1 - 1 / aVals[i]**2) for i in range(len(rVals))]
    derivVals = []

    for i in range(amtVals - 1):
        deriv = (massVals[i + 1] - massVals[i])/(rVals[i+1] - rVals[i])
        derivVals.append(deriv)

    minInd, minVal = localMinima(derivVals)
    if(np.abs(minVal) < 1e-7):
        return True
    return False

# returns the mass at the radius for which M'(r) is at its minimum, the radius where this happens and also the index of that radius
def getGravitationalMass_improved(sol):
    rVals = sol['t']
    aVals = sol['y'][0]
    amtVals = len(rVals)

    massVals = [rVals[i] / 2 * (1 - 1 / aVals[i]**2) for i in range(len(rVals))]
    derivVals = []

    for i in range(amtVals - 1):
        deriv = (massVals[i + 1] - massVals[i])/(rVals[i+1] - rVals[i])
        derivVals.append(deriv)

    minInd, minVal = localMinima(derivVals)
    return massVals[minInd], rVals[minInd], minInd

def newKernel(phic, rhocVals, index):
    rhoc = rhocVals[index]
    omegaRes = binarySearch(phic, rhoc, 1, 2, terminateMass)
    res = getProtoSolution(phic, rhoc, omegaRes)
    mass = getGravitationalMass_improved(res)
    
    tempDict = {"phic" : phic, "rhoc" : rhoc, "sol" : res, "omega" : omegaRes, "mass" : mass, "terminator" : "mass"}
    

    print("finished", phic, rhoc)
    print("mass is", mass)
    return tempDict

def parallelAction(phic, rhocVals):
    print("starting", phic)
    rhocVals = list(rhocVals)

    rhocAmt = len(rhocVals)
    # print("hsdlasd")
    threadPool = multiprocessing.Pool(10)
    # print("hsdlasd")
    result = threadPool.starmap(newKernel, zip(repeat(phic), repeat(rhocVals), range(rhocAmt)))
    # threadPool.map(kernel2, (range(phicamt), ))
    # print("hsdlasd")

    return result

import sys
#import winsound

def newKernel_v2(valTuple):
    phic = valTuple[0]
    rhoc = valTuple[1]
    omegaRes, omegaMin, omegaMax = binarySearch(phic, rhoc, 1, 2, terminateMass)
    res = getProtoSolution(phic, rhoc, omegaRes)
    mass, massRad, massRadIdx = getGravitationalMass_improved(res)
    fNumber = getFermionNumber(res)
    fRadius = getFermionRadius(res)
    bNumber = getBaryonNumber(res, omegaRes, massRadIdx)

    tempDict = {
        "phic" : phic, # code units
        "rhoc" : rhoc, # code units
        # "sol" : res, 
        "omega" : omegaRes, # code units
        "omegaMin": omegaMin, # code units
        "omegaMax": omegaMax, # code units
        "mass" : mass, # sun masses
        "fNumber" : fNumber, # code units
        "bNumber" : bNumber, # code units
        "fRadius" : fRadius, # km
        "Nb/Nf" : bNumber / fNumber, # code units
        "Nf/Nb" : fNumber / bNumber, # code units
        "convergenceRad" : massRad, # code units
        "convergenceIdx" : massRadIdx, # ... 
        "terminator" : "mass"
        }
    
    frequency = 440  # Set Frequency To 2500 Hertz
    duration = 750  # Set Duration To 1000 ms == 1 second
    #winsound.Beep(frequency, duration)

    # sys.stdout.write("finished {} {} mass {}".format(phic, rhoc, mass))
    return tempDict

def parallelAction_v2(allVals):
    threadPool = multiprocessing.Pool(12)
    result = threadPool.map(newKernel_v2, allVals)
    # result = map(newKernel_v2, allVals)
    # threadPool.map(kernel2, (range(phicamt), ))
    # print("hsdlasd")

    return result


################################################
### everything below this line is depricated ###
################################################




# sol, solOmega = getSolution(0.01, 0.05)
# print(solOmega)
def kernel(rhoc, sols, phicVals, index):
    print("asd")
    print("starting", index)
    sys.stdout.flush()

    sol, omega = getSolution(phicVals[index], rhoc)
    # sols[index] = [sol, omega]

    print("finished", index)
    return [sol, omega]

# threadPool.starmap(kernel, zip(repeat(rhoc), repeat(sols), repeat(phicVals), range(phicamt), repeat(phicamt)))


class scalarFieldConfig:
    def __init__(self, sol, omega) -> None:
        self.baryonNumberCount = getBaryonNumber(sol, omega)
        self.fermionNumberCount = getFermionNumber(sol)
        self.totalGravitationalMass = getGravitationalMass(sol)
        self.fermionRadius = getFermionRadius(sol)
        self.omega = omega
        self.sol = sol

    baryonNumberCount = 0
    fermionNumberCount = 0
    totalGravitationalMass = 0
    fermionRadius = 0
    omega = 0
    sol = None


def calc_scalar_parallel(rhoc, phicVals, sols):

    # def kernel2(index):
    #     print("asd")
    #     print("starting", index)
    #     sys.stdout.flush()

    #     sol, omega = getSolution(phicVals[index], rhoc)
    #     sols[index] = [sol, omega]
    #     print("finished", index)


    phicamt = len(phicVals)
    # print("hsdlasd")
    threadPool = multiprocessing.Pool(10)
    # print("hsdlasd")
    result = threadPool.starmap(kernel, zip(repeat(rhoc), repeat(sols), repeat(phicVals), range(phicamt)))
    # threadPool.map(kernel2, (range(phicamt), ))
    # print("hsdlasd")

    finalResult = []

    for entry in result:
        finalResult.append(scalarFieldConfig(entry[0], entry[1]))

    return finalResult



few = 5 * [0]


def test_kernel(index, list):
    print(index)
    few[index] = index
    return index

def run_test_kernel(list):
    
    threadPool = multiprocessing.Pool(1)

    return threadPool.starmap(test_kernel, zip(range(len(list)), repeat(list)))
    # threadPool.map(test_kernel, range(len(list)))

    # x = multiprocessing.Process(target = kernel, args = (rhoc, sols, phicVals, 0, 1, ))
    # x.start()
    # x.join()
# kernel(rhoc, sols, phicVals, 0, 1)


if __name__ == "__main__":



    print(run_test_kernel(few))
    print(few)

    exit(1)

    rhoc = 0.003

    phicmin = 0.01
    phicmax = 0.5
    phicamt = 2

    sols = phicamt * [None]
    phicVals = np.linspace(phicmin, phicmax, phicamt)

    threadPool = multiprocessing.Pool(2)

    threadPool.starmap(kernel, zip(repeat(rhoc), repeat(sols), repeat(phicVals), range(phicamt)))

