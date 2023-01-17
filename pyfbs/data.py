import numpy as np
import os
from .pyfbs_cython import PyFermionBosonStarTLN

## This file is simply intented to convert the output of either the cpp code or the pyfbs into numpy arrays that are used for plotting

# loads the custum FBS data values from a file written by the cpp code
def load_file(filename):
    f = open(filename)
    l0 = next(f)
    labels = l0.replace('#', '').strip().split('\t')
    indices = dict([(labels[i].strip(), i) for i in range(len(labels))])
    f.close()
    data = np.loadtxt(filename) #, delimeter = ' ')
    return data, indices

# Converts a curve of PyFermionBosonStar(TLN) elements into a numpy array
def simplify_curve(PyMR_curve):
    data = []
    indices = {}
    labels = ['M_T', 'rho_0', 'phi_0', 'R_F', 'N_F', 
                'R_B', 'R_B_0', 'N_B', 'omega', 'mu', ]
    if isinstance(PyMR_curve[0], PyFermionBosonStarTLN):
        labels.append('k2')
        labels.append('lambda_tidal')
        labels.append('phi_1_0')
        labels.append('H_0')
    
    for l, i in zip(labels, range(len(labels))):
        indices[l] = i

    for fbs in PyMR_curve:
        g = fbs.get()
        d = []
        for i in range(len(labels)):
            d.append(g[labels[i]])
        data.append(d)
   
    return np.array(data), indices

