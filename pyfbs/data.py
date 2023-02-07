import numpy as np
import os
try:
    from .pyfbs_cython import PyFermionBosonStarTLN, PyFermionBosonStar
except ImportError:
    pass

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

def add_Lambda_int(df, indices):
    if not 'Lambda_int' in indices or indices['Lambda_int'] >= np.shape(df)[1]:
        df = np.append(df, np.transpose([ df[:,indices['lambda']]/df[:,indices['mu']]**2 / 8. / np.pi  ]), axis=1)
        indices['Lambda_int'] = np.shape(df)[1]-1
    return df, indices
    
def add_Lambda_tidal(df, indices):
    if not 'Lambda_tidal' in indices or indices['Lambda_tidal'] >= np.shape(df)[1]:
        df = np.append(df, np.transpose([ df[:,indices['lambda_tidal']]/df[:,indices['M_T']]**5  ]), axis=1)
        indices['Lambda_tidal'] = np.shape(df)[1]-1
    return df, indices
    
def add_R_max(df, indices):
    if not 'R_m' in indices or indices['R_m'] >= np.shape(df)[1]:
        df = np.append(df, np.transpose([ np.maximum(df[:,indices['R_F']], df[:,indices['R_B']] ) ]), axis=1)
        indices['R_m'] = np.shape(df)[1]-1
    return df, indices


# Converts a curve of PyFermionBosonStar(TLN) elements into a numpy array
def simplify_curve(PyMR_curve):
    data = []
    indices = {}
    labels = ['M_T', 'rho_0', 'phi_0', 'R_F', 'N_F', 
                'R_B', 'R_B_0', 'N_B', 'omega', 'mu', 'R_G']
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

def file_to_curve(df, indices, pyEoS):
    PyMR_curve = []
    for d in df:
        PyMR_curve.append(PyFermionBosonStar.FromParameters(pyEoS, d[indices['mu']], lambda_ = d[indices['lambda']], 
                                omega = d[indices['omega']], rho_0 = d[indices['rho_0']], phi_0 = d[indices['phi_0']] ) )            
    return PyMR_curve
    
