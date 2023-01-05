from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr, make_shared
from libcpp cimport bool

from cpython cimport array
import array
cimport numpy as np
import numpy as np


from cpyfbs cimport *

cdef class PyEoS:
    cdef shared_ptr[EquationOfState] eos

    def callEOS(self, P):
        cdef double rho, eps
        deref(self.eos).callEOS(rho, eps, P)
        return rho, eps

    def get_P_from_rho(self, rho_in, epsilon=0.):
        return deref(self.eos).get_P_from_rho(rho_in, epsilon)

    def min_P(self):
        return deref(self.eos).min_P()

    def min_rho(self):
        return deref(self.eos).min_rho()

cdef class PyEoStable(PyEoS):
#    cdef shared_ptr[EoStable] eos

    def __cinit__(self, str filename):
        cdef string f = <string> filename.encode('utf-8')
        cdef shared_ptr[EoStable] ptr =  make_shared[EoStable](f)
        self.eos = ptr
        # see https://stackoverflow.com/questions/67626270/inheritance-and-stdshared-ptr-in-cython


cdef class PyPolytropicEoS(PyEoS):

    def __cinit__(self, double kappa, double Gamma):
        self.eos = make_shared[PolytropicEoS](kappa, Gamma)

cdef class PyCausalEoS(PyEoS):

    def __cinit__(self, double eps_f, double P_f=0.):
        self.eos = make_shared[CausalEoS](eps_f, P_f)


cdef class PyIntegrationOptions:
    cdef shared_ptr[IntegrationOptions] io

    def __cinit__(self, int max_step=1000000, double target_error=1e-12, double min_stepsize=1e-16, double max_stepsize=1e-2, bool save_intermediate=False, int verbose=0):
        self.io = make_shared[IntegrationOptions](max_step, target_error, min_stepsize, max_stepsize, save_intermediate, verbose)



cdef class PyFermionBosonStar:
    cdef shared_ptr[FermionBosonStar] fbs
    cdef bool evaluated
    cdef np.ndarray results

    def __cinit__(self):
        self.evaluated=False

    @staticmethod
    def FromParameters(PyEoS pyEoS, mu, lambda_=0., omega=0., rho_0 = 0., phi_0 = 0.):
        cdef PyFermionBosonStar pfbs = PyFermionBosonStar.__new__(PyFermionBosonStar)
        pfbs.fbs = make_shared[FermionBosonStar](pyEoS.eos, <double>mu, <double>lambda_, <double>omega, <double>rho_0, <double>phi_0)
        return pfbs

    @staticmethod
    cdef PyFermionBosonStar FromPointer(shared_ptr[FermionBosonStar] fbs):
        cdef PyFermionBosonStar pfbs = PyFermionBosonStar.__new__(PyFermionBosonStar)
        pfbs.fbs = fbs
        return pfbs

    @staticmethod
    cdef PyFermionBosonStar FromObject(FermionBosonStar fbs):
        cdef PyFermionBosonStar pfbs = PyFermionBosonStar.__new__(PyFermionBosonStar)
        pfbs.fbs = make_shared[FermionBosonStar](fbs)
        return pfbs

    def bisection(self, omega_0, omega_1, n_mode=0, max_step=500, delta_omega=1e-15):
        deref(self.fbs).bisection(omega_0, omega_1, n_mode, max_step, delta_omega)
        self.evaluated=False


    def evaluate_model(self, PyIntegrationOptions intOpts=PyIntegrationOptions()):
        cdef stdvector[step] res
        cdef string empty
        deref(self.fbs).evaluate_model(res, deref(intOpts.io), empty)
        self.evaluated=True
        cdef unsigned int i,j
        self.results = np.zeros([res.size(), res[0].second.size()+1])
        for i in range(res.size()):
            self.results[i, 0] = res[i].first
            for j in range(res[i].second.size()):
                self.results[i,1+j] = res[i].second[j]
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
        return {"M_T":deref(self.fbs).M_T,
                "N_B":deref(self.fbs).N_B,
                "N_F":deref(self.fbs).N_F,
                "R_B":deref(self.fbs).R_B,
                "R_B_0":deref(self.fbs).R_B_0,
                "R_F":deref(self.fbs).R_F,
                "R_F_0":deref(self.fbs).R_F_0,
                "rho_0":deref(self.fbs).rho_0,
                "phi_0":deref(self.fbs).phi_0,
                "mu":deref(self.fbs).mu,
                "omega":deref(self.fbs).omega,
                }


cdef class PyFermionBosonStarTLN(PyFermionBosonStar):
    cdef shared_ptr[FermionBosonStarTLN] fbstln

    @staticmethod
    def FromFBS(PyFermionBosonStar pfbs_notln):
        cdef PyFermionBosonStarTLN pfbs = PyFermionBosonStarTLN.__new__(PyFermionBosonStarTLN)
        pfbs.fbstln = make_shared[FermionBosonStarTLN](deref(pfbs_notln.fbs))
        pfbs.fbs = pfbs.fbstln
        return pfbs

    @staticmethod
    cdef PyFermionBosonStarTLN FromPointer(shared_ptr[FermionBosonStarTLN] fbstln):
        cdef PyFermionBosonStarTLN pfbs = PyFermionBosonStarTLN.__new__(PyFermionBosonStarTLN)
        pfbs.fbstln = fbstln
        pfbs.fbs = pfbs.fbstln
        return pfbs

    @staticmethod
    cdef PyFermionBosonStarTLN FromObject(FermionBosonStarTLN fbstln):
        cdef PyFermionBosonStarTLN pfbs = PyFermionBosonStarTLN.__new__(PyFermionBosonStarTLN)
        pfbs.fbstln = make_shared[FermionBosonStarTLN](fbstln)
        pfbs.fbs = pfbs.fbstln
        return pfbs

    def set_initial_conditions(self, phi_1_0, H_0):
        #deref(self.fbstln).set_initial_conditions(phi_1_0, H_0)
        deref(self.fbstln).phi_1_0 = phi_1_0
        deref(self.fbstln).H_0 = H_0
        self.evaluated=False

    def bisection_phi_1(self, phi_1_0, phi_1_1, n_mode = 0, max_step = 200, delta_phi_1=1e-10):
        deref(self.fbstln).bisection_phi_1(phi_1_0, phi_1_1, n_mode, max_step, delta_phi_1)
        self.evaluated=False

    def evaluate_model(self):
        cdef stdvector[step] res
        cdef string empty
        deref(self.fbstln).evaluate_model(res, empty)
        self.evaluated=True
        cdef unsigned int i,j
        self.results = np.zeros([res.size(), res[0].second.size()+1])
        for i in range(res.size()):
            self.results[i, 0] = res[i].first
            for j in range(res[i].second.size()):
                self.results[i,1+j] = res[i].second[j]
        return self.results

    def plot(self, ax, components=[0,1,2,3,4,5,6,7,8], label=""):
        component_labels=["a", r"\alpha", r"\phi", r"\Psi", "P", "H", "dH", r"\phi_1", r"d\phi_1"]
        for c in components:
            ax.plot(self.results[:,0], self.results[:,1+c], label=(label + f"${component_labels[c]}$"))

    def get(self):
        d = PyFermionBosonStar.get(self)
        d['k2'] = deref(self.fbstln).k2
        d['lambda_tidal'] = deref(self.fbstln).lambda_tidal
        d['phi_1_0'] = deref(self.fbstln).phi_1_0
        d['H_0'] = deref(self.fbstln).H_0
        d['y_max'] = deref(self.fbstln).y_max
        d['R_ext'] = deref(self.fbstln).R_ext
        return d


cdef class PyMRcurve:

    @staticmethod
    def from_list(pMRphi_curve, str filename=""):
        cdef string f = <string> filename.encode('utf-8')
        cdef stdvector[FermionBosonStar] MRphi_curve
        cdef PyFermionBosonStar fbs

        for fbs in pMRphi_curve:
            MRphi_curve.push_back(deref( fbs.fbs))

        calc_rhophi_curves(MRphi_curve, 2);

        if(not f.empty()):
            write_MRphi_curve(MRphi_curve, f)

        pMRphi_curve = []
        for i in range(MRphi_curve.size()):
            pMRphi_curve.append(PyFermionBosonStar.FromObject(MRphi_curve[i]))

        return pMRphi_curve


    @staticmethod
    def from_rhophi_list(mu, lam, PyEoS eos, np.ndarray rho_c_grid, np.ndarray phi_c_grid, str filename=""):
        cdef stdvector[double] crho_c_grid
        cdef stdvector[double] cphi_c_grid
        cdef stdvector[FermionBosonStar] MRphi_curve
        cdef string f = <string> filename.encode('utf-8')
        cdef int i
        for i in range(len(rho_c_grid)):
            crho_c_grid.push_back(rho_c_grid[i])
        for i in range(len(phi_c_grid)):
            cphi_c_grid.push_back(phi_c_grid[i])

        calc_rhophi_curves(mu, lam, eos.eos, crho_c_grid, cphi_c_grid, MRphi_curve, 2)
        if(not f.empty()):
            write_MRphi_curve(MRphi_curve, f)

        pMRphi_curve = []
        for i in range(MRphi_curve.size()):
                pMRphi_curve.append(PyFermionBosonStar.FromObject(MRphi_curve[i]))
        return pMRphi_curve

    @staticmethod
    def from_NbNf_curve(mu, lam, PyEoS eos, np.ndarray rho_c_grid, np.ndarray NbNf_grid, str filename=""):
        cdef stdvector[double] crho_c_grid
        cdef stdvector[double] cNbNf_grid
        cdef stdvector[FermionBosonStar] MRphi_curve
        cdef string f = <string> filename.encode('utf-8')
        cdef int i
        for i in range(len(rho_c_grid)):
            crho_c_grid.push_back(rho_c_grid[i])
        for i in range(len(NbNf_grid)):
            cNbNf_grid.push_back(NbNf_grid[i])

        calc_NbNf_curves(mu, lam, eos.eos, crho_c_grid, cNbNf_grid, MRphi_curve)

        if(not f.empty()):
            write_MRphi_curve(MRphi_curve, f)

    @staticmethod
    def calc_TLN_curve(pMRphi_curve, str filename=""):
        cdef stdvector[FermionBosonStar] MRphi_curve
        cdef stdvector[FermionBosonStarTLN] tln_curve
        cdef PyFermionBosonStar fbs
        cdef string f = <string> filename.encode('utf-8')

        for fbs in pMRphi_curve:
            MRphi_curve.push_back(deref( fbs.fbs))

        calc_MRphik2_curve(MRphi_curve, tln_curve, 2)

        if(not f.empty()):
            write_MRphi_curve(tln_curve, f)

        pTLN_curve = []
        for i in range(tln_curve.size()):
            pTLN_curve.append(PyFermionBosonStarTLN.FromObject(tln_curve[i]))

        return pTLN_curve

