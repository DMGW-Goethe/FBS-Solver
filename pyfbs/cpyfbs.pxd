
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr


#cdef extern from "vector.cpp":
#    pass
#cdef extern from "vector.hpp":
#    pass


#cdef extern from "eos.cpp":
#    pass

cdef extern from "eos.hpp":
    cdef cppclass EoStable:
        EoStable(const string filename) except +

        void callEOS(double& myrho, double& epsilon, const double P)
        double get_P_from_rho(const double rho_in)


#cdef extern from "nsmodel.cpp":
#    pass

cdef extern from "nsmodel.hpp":
    cdef cppclass FermionBosonStar:
        FermionBosonStar( shared_ptr[EoStable], double, double, double) except +

        double mu
        #double lambda # doesn't work due to variable name
        double omega

        void set_initial_conditions(const double rho_0, const double phi_0)
        void bisection(double omega_0, double omega_1, int n_mode, int max_step, double delta_omega)
        void evaluate_model(string filename)
        void shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, double omega_0, double omega_1, int n_mode, int max_step, double delta_omega)

        double M_T
        double N_B
        double N_F
        double R_B
        double R_F
        double R_F_0



