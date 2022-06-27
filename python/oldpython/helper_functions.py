import numpy as np
import sys
from types import ModuleType, FunctionType
from gc import get_referents

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

# taken from https://stackoverflow.com/a/1751478
def reformat(l, n):
    n = max(1, n)
    return [l[i:i+n] for i in range(0, len(l), n)]

# central (symmetric) stencil
def gradientCentralDifference(vals, i, j, order, prop):

    stencil2nd = [-0.5, 0.0, 0.5]
    stencil4th = [1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0]
    stencil6th = [-1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0]

    stencil = []
    if order == 2:
        stencil = stencil2nd
    elif order == 4:
        stencil = stencil4th
    elif order == 6:
        stencil = stencil6th

    derivIdir = 0
    derivJdir = 0
    for iTemp in range(-order//2, order//2 + 1):
        # print(iTemp)
        derivIdir += prop(vals[i + iTemp][j]) * stencil[iTemp + order//2]
        derivJdir += prop(vals[i][j + iTemp]) * stencil[iTemp + order//2]

    derivIdir /= vals[i + 1][j].phic - vals[i][j].phic
    derivJdir /= vals[i][j + 1].rhoc - vals[i][j].rhoc
    return derivIdir, derivJdir

# forward stenctils
def gradientForwardDifference(vals, i, j, order, prop):
    stencil1st = [-1.0, 1.0]
    stencil2nd = [-3.0/2.0, 2.0, -1.0/2.0]
    stencil3rd = [-11.0/6.0, 3.0, -3.0/2.0, 1.0/3.0]
    stencil4th = [-25.0/12.0, 4.0, -3.0, 4.0/3.0, -1.0/4.0]

    stencil = []
    if order == 1:
        stencil = stencil1st
    elif order == 2:
        stencil = stencil2nd
    elif order == 3:
        stencil = stencil3rd
    elif order == 4:
        stencil = stencil4th

    derivIdir = 0
    derivJdir = 0
    for iTemp in range(0, order + 1):
        # print(iTemp)
        derivIdir += prop(vals[i + iTemp][j]) * stencil[iTemp]
        derivJdir += prop(vals[i][j + iTemp]) * stencil[iTemp]

    derivIdir /= vals[i + 1][j].phic - vals[i][j].phic
    derivJdir /= vals[i][j + 1].rhoc - vals[i][j].rhoc
    return derivIdir, derivJdir


# Following function can calculate the size of custom made objects.
# Taken from https://stackoverflow.com/a/30316760
# Custom objects know their class.
# Function objects seem to know way too much, including modules.
# Exclude modules as well.
BLACKLIST = type, ModuleType, FunctionType

def getsize(obj):
    """sum size of object & members."""
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: '+ str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size
