import numpy as np
from numba import jit

import helper_functions
import integrator

import importlib
importlib.reload(helper_functions)
importlib.reload(integrator)

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

    minSol = integrator.getProtoSolution(phic, rhoc, omegaMin)
    maxSol = integrator.getProtoSolution(phic, rhoc, omegaMax)

    minSolPhi = minSol['y'][2]
    maxSolPhi = maxSol['y'][2]

    # minSol should not cross phi = 0 at any point
    # maxSol should only cross once, so as a initial check we can decrease maxSol until only one intersection with the r axis is left
    if(helper_functions.amtOfRoots(minSolPhi) != 0):
        print("Error: the lower bound for omega should result in a distribution of phi with no roots at all")
        return 0.0, 0.0, 0.0

    # if the upper value of omega does not result in a root that is being sandwhiched, then we increase it until there is one
    while(helper_functions.amtOfRoots(maxSolPhi) == 0):
        omegaMin = omegaMax
        omegaMax += 0.1

        minSol = integrator.getProtoSolution(phic, rhoc, omegaMin)
        maxSol = integrator.getProtoSolution(phic, rhoc, omegaMax)

        minSolPhi = minSol['y'][2]
        maxSolPhi = maxSol['y'][2]

    # if(amtOfRoots(maxSolPhi) == 0):
    #     print("Error: the upper bound of omega should result in at least one root")
    #     return 0
    
    
    # Part1: decrease omegaMax until only the ground state is sandwiched by omegaMin and omegaMax
    # for this I define two temporary variables that will keep track on the upper and lower limit for omegaMax
    omegaMaxL = omegaMin
    omegaMaxU = omegaMax

    currAmtRoots = helper_functions.amtOfRoots(maxSolPhi) # this is how may roots are currently present for omega = omegaMax

    while True:
        # if we currently already have only one root, then we do not have to adjust omegaMax and can just skip this part
        if(currAmtRoots == 1):
            break
        
        # now we will test the value sad is in the middle of omegaMaxL and omegaMaxU and see how many roots are enclosed
        omegaMaxTemp = (omegaMaxL + omegaMaxU) / 2.0
        maxSolTemp = integrator.getProtoSolution(phic, rhoc, omegaMaxTemp)
        maxSolPhiTemp = maxSolTemp['y'][2]
        
        if(helper_functions.amtOfRoots(maxSolPhiTemp) == 0): # if we have zero roots, then the guess for omegaMax was too low and we have to increase it again. For this we can just set omegaMaxL as this value and continue
            omegaMaxL = omegaMaxTemp
            continue
        elif(helper_functions.amtOfRoots(maxSolPhiTemp) == 1): # if we have exactly one root then we are done
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

    minSol = integrator.getProtoSolution(phic, rhoc, omegaMin)
    maxSol = integrator.getProtoSolution(phic, rhoc, omegaMax)

    for i in range(maxIteration):
        # the new guess will be generated from the value in the middle of omegaMin and omegaMax
        omegaMid = (omegaMin + omegaMax) / 2.0
        midSol = integrator.getProtoSolution(phic, rhoc, omegaMid)
        
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

    minInd, minVal = helper_functions.localMinima(derivVals)
    if(np.abs(minVal) < 1e-7):
        return True
    return False