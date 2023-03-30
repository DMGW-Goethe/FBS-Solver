import numpy as np
import matplotlib.pyplot as plt
import pyfbs_cython as pyfbs


eosDD2 = pyfbs.PyEoStable("../EOS_tables/eos_HS_DD2_with_electrons.beta")

rho_cs = 5e-3 * np.linspace(0., 1., 120)**3     # in m_planck
phi_cs = 0.09 * np.linspace(0., 1., 120)**3     # in m_planck
mus = np.array([.1, 1., 10.])                   # in m_planck
lams = np.array([0., 10., 100.])                # Lambda_int

filedir = "../data/"

for mu in mus:
    for lam in lams:
        pMR = []
        filename = filedir + "mu" + str(mu) + "_lam" + str(lam) + "_" + str(len(rho_cs)) + "x" + str(len(phi_cs)) + ".txt"
        for rho_c in rho_cs:
            for phi_c in phi_cs:
                pMR.append(pyfbs.PyFermionBosonStar.FromParameters(eosDD2, mu, lambda_=lam*8.*np.pi*mu**2, rho_0=rho_c, phi_0=phi_c))
        print(filename)
        pMR = pyfbs.PyMRcurve.from_list(pMR, filename)
        tln_curve = pyfbs.PyMRcurve.calc_TLN_curve(pMR, filename)


