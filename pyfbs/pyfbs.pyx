from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr, make_shared
from libcpp cimport bool

cimport numpy as np
import numpy as np


from cpyfbs cimport *


cdef class PyEoStable:
    cdef shared_ptr[EoStable] eos

    def __cinit__(self, str filename):
        cdef string f = <string> filename.encode('utf-8')
        self.eos = make_shared[EoStable](f)

#    def __dealloc__(self):
#        del self.eos

    def callEOS(self, P):
        cdef double rho, eps
        deref(self.eos).callEOS(rho, eps, P)
        return rho, eps

    def get_P_from_rho(self, rho_in):
        return deref(self.eos).get_P_from_rho(rho_in)

cdef class PyPolytropicEoS:
    cdef shared_ptr[PolytropicEoS] eos

    def __cinit__(self, double kappa, double Gamma):
        self.eos = make_shared[PolytropicEoS](kappa, Gamma)

#    def __dealloc__(self):
#        del self.eos

    def callEOS(self, P):
        cdef double rho, eps
        deref(self.eos).callEOS(rho, eps, P)
        return rho, eps

    def get_P_from_rho(self, rho_in):
        return deref(self.eos).get_P_from_rho(rho_in)


cdef class PyIntegrationOptions:
    cdef shared_ptr[IntegrationOptions] io

    def __cinit__(self, int max_step=1000000, double target_error=1e-10, double min_stepsize=1e-18, double max_stepsize=1e-2, bool save_intermediate=False, int verbose=0):
        self.io = make_shared[IntegrationOptions](max_step, target_error, min_stepsize, max_stepsize, save_intermediate, verbose)



cdef class PyFermionBosonStar:
    cdef shared_ptr[FermionBosonStar] fbs
    cdef bool evaluated
    cdef np.ndarray results

    def __cinit__(self, PyEoStable pyEoS, mu, lambda_=0., omega=0.):
        self.fbs = make_shared[FermionBosonStar](pyEoS.eos, <double>mu, <double>lambda_, <double>omega)
        self.evaluated=False

    def set_initial_conditions(self, rho_0, phi_0):
        deref(self.fbs).set_initial_conditions(rho_0, phi_0)
        self.evaluated=False

    def bisection(self, omega_0, omega_1, n_mode=0, max_step=500, delta_omega=1e-15):
        deref(self.fbs).bisection(omega_0, omega_1, n_mode, max_step, delta_omega)
        self.evaluated=False

#    def evaluate_model(self):
#        deref(self.fbs).evaluate_model()
#        self.evaluated=True

    def evaluate_model(self):
        cdef stdvector[step] res
        deref(self.fbs).evaluate_model(res)
        self.evaluated=True
        cdef unsigned int i,j
        self.results = np.zeros([res.size(), res[0].second.size()+1])
        for i in range(res.size()):
            self.results[i, 0] = res[i].first
            for j in range(res[i].second.size()):
                self.results[i,j] = res[i].second[j]
        return self.results

    def shooting_NbNf_ratio(self, NbNf_ratio, NbNf_accuracy, omega_0, omega_1, n_mode=0, max_step=500, delta_omega=1e-15):
        deref(self.fbs).shooting_NbNf_ratio(NbNf_ratio, NbNf_accuracy, omega_0, omega_1, n_mode, max_step, delta_omega)
        self.evaluated=True

    def makeFBSstar(self, rho_0, NbNf_ratio, omega_0=1., omega_1=20., r_init=1e-6, n_mode=0, max_step=500, delta_omega=1e-15, NbNf_accuracy=1e-4):
        self.set_initial_conditions(rho_0, 0.)
        self.shooting_NbNf_ratio(NbNf_ratio, NbNf_accuracy, omega_0=omega_0, omega_1=omega_1, n_mode=n_mode, max_step=max_step, delta_omega=delta_omega)
        self.evaluate_model()

    def plot(self, ax, components=[0,1,2,3,4], label=""):
        component_labels=["a", r"\alpha", r"\Phi", r"\Psi", "P"]
        for c in components:
            ax.plot(self.results[:,0], self.results[:,1+c], label=(label + f"${component_labels[c]}$"))


    def get(self):
        if not self.evaluated:
            return
        return {"M_T":deref(self.fbs).M_T,
                "N_B":deref(self.fbs).N_B,
                "N_F":deref(self.fbs).N_F,
                "R_B":deref(self.fbs).R_B,
                "R_F":deref(self.fbs).R_F,
                "R_F_0":deref(self.fbs).R_F_0,
                "rho_0":deref(self.fbs).rho_0,
                "phi_0":deref(self.fbs).phi_0,
                }

