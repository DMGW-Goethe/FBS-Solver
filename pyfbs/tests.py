import numpy as np
import matplotlib.pyplot as plt
from . import data
from .pyfbs_cython import PyFermionBosonStar, PyMRcurve

def generate_testset(filename):
    rho_cs = np.linspace(0., 1e-2, 5)**3
    phi_cs = np.linspace(0.,  1e-1, 5)**3
    mus = np.geomspace(1e-1, 1e1, 3)
    lam_ints = np.append([0.], np.geomspace(1e1, 1e3, 3))

    pMR = []
    for rho_c in rho_cs:
        for phi_c in phi_cs:
            for mu in mus:
                for lam_int in lam_ints:
                    pMR.append(PyFermionBosonStar.FromParameters(eosDD2, mu, lambda_=lam_int*8*np.pi*mu**2, rho_0=rho_c, phi_0=phi_c))
    pMR = PyMRcurve.from_list(pMR)
    tln_curve = PyMRcurve.calc_TLN_curve(pMR, filename)
    return tln_curve


def plot_differences(ax_val, ax_rel, data, indices, ref_data, ref_indices, toplot, plot_ref=False, label=""):
    ref_val = ref_data[:, ref_indices[toplot]]
    #ref_val = ref_data[np.where(ref_data[:, indices['phi_0']] > 1e-5)[0], indices[toplot]]

    if plot_ref:
        ax_val.plot(ref_val, 'o', label='ref')

    val = data[:, indices[toplot]]
    #print(np.where(ref_data[:, indices['phi_0']] > 1e-5)[0])
    #val =         data[np.where(ref_data[:, indices['phi_0']] > 1e-5)[0], indices[toplot]]
    l, = ax_val.plot(val, 'o', label=label)
    rel_err = np.abs(val-ref_val)/ref_val
    ax_rel.plot(rel_err, color=l.get_c())
    rel_err= np.nan_to_num(rel_err)
    print(toplot, np.max(rel_err), data[np.argmax(rel_err), :], label)


def plot(files, toplot, ref_data, ref_indices):
    fig, ax = plt.subplots(2,1)

    for f,i in zip(files, range(len(files))):
        df, indices = data.load_file(f + '.txt')
        plot_differences(ax[0], ax[1], df, indices, ref_data, ref_indices, toplot, plot_ref=(i == 0), label=f)
    ax[0].legend()
    ax[0].grid(); ax[1].grid()
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    #ax[1].set_ylim(-0.001, 0.001)
    ax[0].set_title(toplot)




def main():
    d = "data/"
    ref_file = d + 'targerr1e-15_maxstep5e-3'
    files =  [ d + 'targerr1e-14_maxstep1e-2',]# 'maxR500_targerr1e-14_maxstep5e-3', 'maxR500_targerr1e-15_maxstep1e-2',]
    toplots = ['M_T', 'lambda_tidal', 'R_F', 'R_B_0',]# 'rho_0', 'phi_0', 'lambda', 'mu']

    ref_data, ref_indices = data.load_file(ref_file + '.txt')

    for toplot in toplots:
        plot(files, toplot, ref_data, ref_indices)

    plt.show()

if __name__ == '__main__':
    main()



