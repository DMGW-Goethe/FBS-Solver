from cProfile import label
import pickle
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from scipy.spatial import ConvexHull
from pathos.multiprocessing import ProcessingPool as Pool

import solution_class
import eigenvalue_finder
import integrator
import stabilityCurve
import helper_functions
import eos_importer
import parameters

import importlib
importlib.reload(solution_class)
importlib.reload(eigenvalue_finder)
importlib.reload(integrator)
importlib.reload(stabilityCurve)
importlib.reload(helper_functions)
importlib.reload(parameters)


def plotStabilityCurve():
    with open("DD2_40x40", "rb") as fp:   # Unpickling
        allSols = pickle.load(fp)


    for sol in allSols:
        sol.EoS = parameters.EoS_v

    rhoToSaturationDensity = 2680 # TODO: I should double check this conversion factor

    phicValsAmt = 40
    rhocValsAmt = 40

    phicmin = 0.0
    phicmax = 0.14
    phicVals = np.linspace(phicmin, phicmax, phicValsAmt)

    rhocmin = 0.0 #* rhoToSaturationDensity
    rhocmax = 0.004# * rhoToSaturationDensity
    rhocVals = np.linspace(rhocmin, rhocmax, rhocValsAmt)

    # print(allSols[22*22].getFermionNumber())
    # print(allSols[22*22].getFermionRadius())

    # print(allSols[12*11].getFermionRadius())
    # radVals = allSols[12*11].radVals()
    # pVals = allSols[12*11].pVals()

    # plt.plot(radVals, pVals)
    # plt.yscale('log')
    # plt.show()

    # exit()

    # fnVals = [res.getFermionRadius() for res in allSols]
    # fnVals = np.array_split(fnVals, phicValsAmt)

    # print(fnVals)
    
    # X, Y = np.meshgrid(rhocVals, phicVals)

    # plt.figure(dpi=120)
    # heatmap = plt.imshow(fnVals, extent=[rhocmin, rhocmax, phicmin, phicmax], origin='lower', cmap='turbo', aspect='auto', interpolation='none', alpha=0.8)
    # # plt.clim(0.0,2.4)
    # cbar = plt.colorbar()
    # cbar.set_label(r"""Total Gravitational Mass  [M$_\odot$]""", rotation=270)
    # cbar.ax.get_yaxis().labelpad = 15
    # plt.show()

    # exit()


    massVals = [res.getGravitationalMass() for res in allSols]
    massVals = np.array_split(massVals, phicValsAmt)


    X, Y = np.meshgrid(rhocVals, phicVals)

    plt.figure(dpi=120)
    heatmap = plt.imshow(massVals, extent=[rhocmin, rhocmax, phicmin, phicmax], origin='lower', cmap='turbo', aspect='auto', interpolation='none', alpha=0.8)
    # plt.clim(0.0,2.4)
    cbar = plt.colorbar()
    cbar.set_label(r"""Total Gravitational Mass  [M$_\odot$]""", rotation=270)
    cbar.ax.get_yaxis().labelpad = 15

    contours = plt.contour(X, Y, massVals, '-',colors=['purple', 'brown', 'red', 'green', 'orange'], levels=[0.62, 1.2, 1.6, 1.9, 2.3])
    for c in contours.collections:
        c.set_dashes([(0, (2.0, 2.0))])
    plt.clabel(contours, inline=True, fontsize=8)
    
    stencilOrder = 3
    massDerivsStuff, fNstuff = stabilityCurve.calcCurve(helper_functions.reformat(allSols, phicValsAmt), stencilOrder)

    X, Y = np.meshgrid(rhocVals[0:-stencilOrder], phicVals[0:-stencilOrder])
    contours = plt.contour(X, Y, massDerivsStuff, colors='black', levels=[0.0])
    plt.clabel(contours, inline=False, fontsize=0)
    plt.xlabel(r"""$\rho_c$ [Saturation Density]""")
    plt.ylabel(r"""$\varphi_c$ [Code Units]""")

    plt.title("DD2 EoS")
    # xtick_temp = plt.xticks()
    # print(xtick_temp)
    # plt.xticks(xtick_temp[0], [tick * rhoToSaturationDensity for tick in xtick_temp[0]])
    p = contours.collections[0].get_paths()
    with open("DD2_40x40_stability_contours", "wb") as fp:   #Pickling
        pickle.dump(p, fp)

    plt.show()



    # # ---- #
    # X, Y = np.meshgrid(rhocVals[0:-3], phicVals[0:-3])


    # plt.figure(dpi=120)
    # plt.imshow(massDerivsStuff, extent=[rhocmin, rhocmax, phicmin, phicmax], origin='lower', cmap='turbo', aspect='auto', interpolation='none')
    # plt.clim(-10, 10)
    # plt.colorbar()

    # contours = plt.contour(X, Y, massDerivsStuff, colors='black', levels=[-1.0, 0.0, 1.0])
    # plt.clabel(contours, inline=True, fontsize=8)
    # plt.show()
    # for kek in b:
        # print(kek.getGravitationalMass())

def scatterPlotter():
    # plotStabilityCurve()

    with open("DD2_40x40_stability_contours", "rb") as fp:   # Unpickling
        contours = pickle.load(fp)

    # for kek in contours:
    #     print(len(kek.vertices))

    # exit()
    yes = contours[3].vertices
    x = yes[:,0]
    y = yes[:,1]

    llim = 0
    ulim = -1

    x = x[llim:ulim]
    y = y[llim:ulim]

    yes = yes[llim:ulim]

    epsilon = 1e-5
    yes = np.insert(yes, 0, [-epsilon, yes[0][1] + epsilon], axis=0)
    yes = np.insert(yes, 0, [-epsilon, -epsilon], axis=0)

    yes = np.append(yes, [[yes[-1][0] + epsilon, -epsilon]], axis=0)

    # print(yes)

    poly = Polygon(yes)
    # testPoint = Point(0.0, 0.0)
    # print(testPoint.within(poly))

    # plt.plot(x,y)
    # plt.show()

    # exit()

    with open("DD2_40x40", "rb") as fp:   # Unpickling
        allSols = pickle.load(fp)
    for sol in allSols:
        sol.EoS = parameters.EoS_v
    # phicValsAmt = 50
    # rhocValsAmt = 50

    # phicmin = 0.0
    # phicmax = 0.14
    # phicVals = np.linspace(phicmin, phicmax, phicValsAmt)

    # rhocmin = 0.0 #* rhoToSaturationDensity
    # rhocmax = 0.008# * rhoToSaturationDensity
    # rhocVals = np.linspace(rhocmin, rhocmax, rhocValsAmt)

    # fnVals = [res.getFermionRadius() for res in allSols]
    # fnVals = np.array_split(fnVals, phicValsAmt)

    # print(fnVals)
    
    # X, Y = np.meshgrid(rhocVals, phicVals)

    # plt.figure(dpi=120)
    # heatmap = plt.imshow(fnVals, extent=[rhocmin, rhocmax, phicmin, phicmax], origin='lower', cmap='turbo', aspect='auto', interpolation='none', alpha=0.8)
    # # plt.clim(0.0,2.4)
    # cbar = plt.colorbar()
    # cbar.set_label(r"""Total Gravitational Mass  [M$_\odot$]""", rotation=270)
    # cbar.ax.get_yaxis().labelpad = 15
    # plt.show()

    # exit()

    allMasses = []
    allRadii = []
    allRadiiAndMass = []

    for sol in allSols:
        if(sol.rhoc == 0):
            continue
        testPoint = Point(sol.rhoc, sol.phic)

        if(testPoint.within(poly)):
            if(sol.getFermionRadius() == -1):
                continue
            allMasses.append(sol.getGravitationalMass())
            allRadii.append(sol.getFermionRadius())

        allRadiiAndMass.append([sol.getFermionRadius(), sol.getGravitationalMass()])
    
    plt.scatter(allRadii, allMasses, alpha=0.4, label="Stable Configurations")


    phicValsAmt = 1
    rhocValsAmt = 100

    phicmin = 0.0
    phicmax = 0.0
    phicVals = np.linspace(phicmin, phicmax, phicValsAmt)

    rhocmin = 0.0002
    rhocmax = 0.0022
    rhocVals = np.linspace(rhocmin, rhocmax, rhocValsAmt)

    allVals = [(phic, rhoc) for phic in phicVals for rhoc in rhocVals]

    # testVals = [[0.04, 0.002], [0.04, 0.001]]

    # print(allVals)
    # print(reformat(allVals, phicValsAmt))
    results = parallelAction_v2(allVals, 12)

    # results = results[10:-1]

    allMasses = []
    allRadii = []
    for sol in results:
        allMasses.append(sol.getGravitationalMass())
        allRadii.append(sol.getFermionRadius())

    
    plt.plot(allRadii, allMasses, 'k-', lw=3, label="Pure NS")


    plt.xlim([0.0, 20.0])
    plt.ylim([0.0, 3.0])
    plt.xlabel(r"""Radius [km]""")
    plt.ylabel(r"""Mass [M$_\odot$]""")
    plt.title("DD2 EoS")
    plt.legend()
    # hull = ConvexHull(allRadiiAndMass)
    # # for simplex in hull.simplices:
    # #     plt.plot(allRadiiAndMass[simplex][0], allRadiiAndMass[simplex][1], 'k-')
    # # plt.plot(allRadiiAndMass[hull.vertices,0], allRadiiAndMass[hull.vertices][1], 'r--', lw=2)
    # # plt.plot(allRadiiAndMass[hull.vertices[0]][0], allRadiiAndMass[hull.vertices[0]][1], 'ro')


    # point1 = allRadiiAndMass[hull.vertices[0]]
    # point2 = allRadiiAndMass[hull.vertices[1]]

    # tpoint1 = [point1[0], point2[0]]
    # tpoint2 = [point1[1], point2[1]]
    # plt.plot(tpoint1, tpoint2, 'ro')
    #     # print(hull.vertices)

    plt.show()

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

    

    # plotStabilityCurve()
    scatterPlotter()
    exit()




    plt.figure(dpi=120)

    # plt.plot([0.0, 10.0], [0.0, 2.0], 'r-', label=r"""$\phi_c = 0.07$""")
    # plt.xlabel(r"""$r$ [km]""")
    # plt.ylabel(r"""Energy Density [Code Units]""")
    # plt.xlim([0, 15])
    # plt.legend()
    # plt.show()
    # exit()

    currConfiguration = solution_class.FieldConfiguration(0.0, 0.002)
    currConfiguration.findSolution()

    rVals = currConfiguration.radVals()
    rVals = 1.477 * rVals

    eDensity = parameters.EoS_v.get_totalEnergyDensity(currConfiguration.pVals())

    plt.plot(rVals, eDensity, 'k-', label=r"""$\phi_c = 0.0$""")

    currConfiguration = solution_class.FieldConfiguration(0.045, 0.002)
    currConfiguration.findSolution()

    rVals = currConfiguration.radVals()
    rVals = 1.477 * rVals

    eDensity = parameters.EoS_v.get_totalEnergyDensity(currConfiguration.pVals())

    plt.plot(rVals, eDensity, 'r-', label=r"""$\phi_c = 0.045$""")

    bDensity = [0.5 * phi**2 for phi in currConfiguration.phiVals()]
    plt.plot(rVals, bDensity, 'r--')

    
    currConfiguration = solution_class.FieldConfiguration(0.09, 0.002)
    currConfiguration.findSolution()

    rVals = currConfiguration.radVals()
    rVals = 1.477 * rVals
    eDensity = parameters.EoS_v.get_totalEnergyDensity(currConfiguration.pVals())

    plt.plot(rVals, eDensity, 'b-', label=r"""$\phi_c = 0.09$""")

    bDensity = [0.5 * phi**2 for phi in currConfiguration.phiVals()]
    plt.plot(rVals, bDensity, 'b--')

    plt.xlabel(r"""$r$ [km]""")
    plt.ylabel(r"""Energy Density [Code Units]""")
    plt.xlim([0, 15])
    plt.legend()
    plt.show()