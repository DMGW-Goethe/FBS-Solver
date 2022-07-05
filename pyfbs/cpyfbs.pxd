from libcpp cimport bool
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libcpp.utility cimport pair
from libcpp.vector cimport vector as stdvector


cdef extern from "vector.hpp":
    cdef cppclass vector:
        ctypedef size_t size_type
        double& operator[](size_type)
        size_type size()

cdef extern from "eos.hpp":
    cdef cppclass EoStable:
        EoStable(const string filename) except +

        void callEOS(double& myrho, double& epsilon, const double P)
        double get_P_from_rho(const double rho_in)

    cdef cppclass PolytropicEoS:
        PolytropicEoS(const double kappa, const double Gamma)

        void callEOS(double& myrho, double& epsilon, const double P)
        double get_P_from_rho(const double rho_in)


cdef extern from "integrator.hpp" namespace "integrator":
    ctypedef pair[double, vector] step
    ctypedef bool (*event_condition)(const double r, const double dr, const vector& y, const vector& dy, const void*params)

    cdef cppclass Event:
        event_condition condition
        bool stopping_condition
        stdvector[step] steps
        bool active
        string name
        Event(event_condition condition, bool stopping_condition, string name)
        void reset()

    cdef cppclass IntegrationOptions:
        int max_step
        double target_error
        double min_stepsize
        double max_stepsize
        bool save_intermediate
        int verbose
        IntegrationOptions(const int max_step, const double target_error, const double min_stepsize, const double max_stepsize, const bool save_intermediate, const int verbose)



cdef extern from "nsmodel.hpp":
    cdef cppclass FermionBosonStar:
        FermionBosonStar( shared_ptr[EoStable], double, double, double) except +

        double mu
        #double lambda # doesn't work due to variable name
        double omega
        double rho_0
        double phi_0

        void set_initial_conditions(const double rho_0, const double phi_0)
        void bisection(double omega_0, double omega_1, int n_mode, int max_step, double delta_omega)
        int integrate(stdvector[step]& result, stdvector[Event]& events, IntegrationOptions intOpts, double r_init, double r_end)
        void evaluate_model()
        void evaluate_model(stdvector[step]& results)
        void shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, double omega_0, double omega_1, int n_mode, int max_step, double delta_omega)

        double M_T
        double N_B
        double N_F
        double R_B
        double R_F
        double R_F_0

        const Event M_converged
        const Event Psi_diverging
        const Event phi_negative
        const Event phi_positive



