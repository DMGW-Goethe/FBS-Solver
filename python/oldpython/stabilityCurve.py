import numpy as np

import helper_functions

import importlib
importlib.reload(helper_functions)

def calcCurve(allSols, stencilOrder):

    phicValsAmt = len(allSols)
    rhocValsAmt = len(allSols[0])

    # print(allSols[10][10].getGravitationalMass(), allSols[10][10].rhoc)
    # print(allSols[10][11].phic, allSols[10][11].rhoc)
    # exit()
    massDerivs4 = []

    # getFermionNumber = scalar_multiproc.getFermionNumber

    fermNumbers = []

    # stencilOrder = 3

    for iPhi in range(0, phicValsAmt - stencilOrder):
        tempVals = []
        tempFn = []
        for iRho in range(0, rhocValsAmt - stencilOrder):
            gradM = [0, 0]
            # gradM[0], gradM[1] = gradientCentralDifference(allSols, iPhi, iRho, stencilOrder, "mass")
            gradM[0], gradM[1] = helper_functions.gradientForwardDifference(allSols, iPhi, iRho, stencilOrder, lambda x : x.getGravitationalMass())

            perp = np.array([gradM[1], -gradM[0]])
            perp /= np.linalg.norm(perp)

            gradFn = [0, 0]
            # gradFn[0] = (getFermionNumber(solP['sol']) - getFermionNumber(solPm['sol']))/(solP['phic'] - solPm['phic'])
            # gradFn[1] = (getFermionNumber(solR['sol']) - getFermionNumber(solRm['sol']))/(solR['rhoc'] - solRm['rhoc'])

            # gradFn[0], gradFn[1] = gradientCentralDifference(allSols, iPhi, iRho, stencilOrder, 'fNumber')
            gradFn[0], gradFn[1] = helper_functions.gradientForwardDifference(allSols, iPhi, iRho, stencilOrder, lambda x : x.getFermionNumber())


            tempVals.append(perp.dot(gradFn))
            # tempFn.append(getFermionNumber(sol0['sol']))
        massDerivs4.append(tempVals)
        fermNumbers.append(tempFn)
    
    return massDerivs4, fermNumbers

def calcCurve2(allSols, stencilOrder):

    phicValsAmt = len(allSols)
    rhocValsAmt = len(allSols[0])

    # print(allSols[10][10].getGravitationalMass(), allSols[10][10].rhoc)
    # print(allSols[10][11].phic, allSols[10][11].rhoc)
    # exit()
    massDerivs4 = []

    # getFermionNumber = scalar_multiproc.getFermionNumber

    fermNumbers = []

    # stencilOrder = 3

    for iPhi in range(stencilOrder, phicValsAmt - stencilOrder):
        tempVals = []
        tempFn = []
        for iRho in range(stencilOrder, rhocValsAmt - stencilOrder):
            gradM = [0, 0]
            # gradM[0], gradM[1] = gradientCentralDifference(allSols, iPhi, iRho, stencilOrder, "mass")
            gradM[0], gradM[1] = helper_functions.gradientCentralDifference(allSols, iPhi, iRho, stencilOrder, lambda x : x.getGravitationalMass())

            perp = np.array([gradM[1], -gradM[0]])
            perp /= np.linalg.norm(perp)

            gradFn = [0, 0]
            # gradFn[0] = (getFermionNumber(solP['sol']) - getFermionNumber(solPm['sol']))/(solP['phic'] - solPm['phic'])
            # gradFn[1] = (getFermionNumber(solR['sol']) - getFermionNumber(solRm['sol']))/(solR['rhoc'] - solRm['rhoc'])

            # gradFn[0], gradFn[1] = gradientCentralDifference(allSols, iPhi, iRho, stencilOrder, 'fNumber')
            gradFn[0], gradFn[1] = helper_functions.gradientCentralDifference(allSols, iPhi, iRho, stencilOrder, lambda x : x.getFermionNumber())


            tempVals.append(perp.dot(gradFn))
            # tempFn.append(getFermionNumber(sol0['sol']))
        massDerivs4.append(tempVals)
        fermNumbers.append(tempFn)
    
    return massDerivs4, fermNumbers

    
