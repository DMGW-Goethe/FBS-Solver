# import multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
import numpy as np
import pickle

import solution_class
import eigenvalue_finder
import integrator
import importlib
importlib.reload(solution_class)
importlib.reload(eigenvalue_finder)
importlib.reload(integrator)

def Kernel(valTuple):
    phic = valTuple[0]
    rhoc = valTuple[1]

    currConfiguration = solution_class.FieldConfiguration(phic, rhoc)
    currConfiguration.findSolution()
    # currConfiguration.sol = None # we might want to delete the profiles and only save the relevant values like omega and such

    print(valTuple)
    return currConfiguration

def parallelAction_v2(allVals, amtThreads):
    threadPool = Pool(amtThreads)
    result = threadPool.map(Kernel, allVals)

    return result

if __name__ == "__main__":

    phicValsAmt = 40
    rhocValsAmt = 40

    phicmin = 0.0
    phicmax = 0.14
    phicVals = np.linspace(phicmin, phicmax, phicValsAmt)

    rhocmin = 0.0
    rhocmax = 0.004
    rhocVals = np.linspace(rhocmin, rhocmax, rhocValsAmt)

    allVals = [(phic, rhoc) for phic in phicVals for rhoc in rhocVals]

    # testVals = [[0.04, 0.002], [0.04, 0.001]]

    # print(allVals)
    # print(reformat(allVals, phicValsAmt))
    results = parallelAction_v2(allVals, 12)

    for res in results:
        print(res.getGravitationalMass())
        res.EoS = None

    with open("DD2_41x40", "wb") as fp:   #Pickling
        pickle.dump(results, fp)