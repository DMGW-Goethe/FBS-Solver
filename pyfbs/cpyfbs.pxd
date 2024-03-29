from libcpp cimport bool
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libcpp.utility cimport pair
from libcpp.vector cimport vector as stdvector


cdef extern from "vector.hpp" namespace "FBS":
    cdef cppclass vector:
        ctypedef size_t size_type
        double& operator[](size_type)
        size_type size()

cdef extern from "eos.hpp" namespace "FBS":
    cdef cppclass EquationOfState:
        void callEOS(double& myrho, double& epsilon, const double P)
        double get_P_from_rho(const double rho_in, const double epsilon)
        double min_P()
        double min_rho()

    cdef cppclass EoStable(EquationOfState):
        EoStable(const string filename) except +

    cdef cppclass PolytropicEoS(EquationOfState):
        PolytropicEoS(const double kappa, const double Gamma)

    cdef cppclass CausalEoS(EquationOfState):
        CausalEoS(const double eps_f, const double P_f)


cdef extern from "integrator.hpp" namespace "FBS::integrator":
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



cdef extern from "nsmodel.hpp" namespace "FBS":
    cdef cppclass FermionBosonStar:
        FermionBosonStar( shared_ptr[EoStable], double mu, double lam, double omega, double rho_0, double phi_0) except +

        double mu
        #double lambda # doesn't work due to variable name
        double omega
        double rho_0
        double phi_0

        void get_initial_conditions()
        int bisection(double omega_0, double omega_1, int n_mode, int max_step, double delta_omega)
        int integrate(stdvector[step]& result, stdvector[Event]& events, IntegrationOptions intOpts, double r_init, double r_end)
        void evaluate_model()
        void evaluate_model(stdvector[step]& results, IntegrationOptions intOpts, string filename)
        void shooting_NbNf_ratio(double NbNf_ratio, double NbNf_accuracy, double omega_0, double omega_1, int n_mode, int max_step, double delta_omega)

        double M_T
        double N_B
        double N_F
        double R_B
        double R_B_0
        double R_F
        double R_F_0
        double R_G

        const Event M_converged
        const Event Psi_diverging
        const Event phi_negative
        const Event phi_positive

    cdef cppclass FermionBosonStarTLN(FermionBosonStar):
        #FermionBosonStarTLN( shared_ptr[EoStable], double, double, double) except +
        FermionBosonStarTLN(const FermionBosonStar& fbs) except +

        double H_0
        double phi_1_0
        double k2
        double lambda_tidal
        double y_max
        double R_ext

        void get_initial_conditions(const double r_init)
        void evaluate_model(stdvector[step]& results, string filename)
        void evaluate_model()
        int bisection_phi_1(double phi_1_0, double phi_1_1, int n_mode, int max_step, double delta_phi_1)

        const Event dphi_1_diverging
        const Event phi_1_negative
        const Event phi_1_positive


cdef extern from "mr_curves.hpp" namespace "FBS":
    void write_MRphi_curve(const stdvector[FermionBosonStar]& MRphi_curve, string filename);
    void write_MRphi_curve(const stdvector[FermionBosonStarTLN]& MRphi_curve, string filename);

    void calc_rhophi_curves(stdvector[FermionBosonStar]& MRphi_curve, int verbose);

    void calc_rhophi_curves(double mu, double lam, shared_ptr[EquationOfState] EOS, const stdvector[double]& rho_c_grid, const stdvector[double]& phi_c_grid, stdvector[FermionBosonStar]& MRphi_curve, int verbose);

    void calc_NbNf_curves(double mu, double lam, shared_ptr[EquationOfState] EOS, const stdvector[double]& rho_c_grid, const stdvector[double]& NbNf_grid, stdvector[FermionBosonStar]& MRphi_curve);

    void calc_MRphik2_curve(const stdvector[FermionBosonStar]& MRphi_curve,  stdvector[FermionBosonStarTLN]& MRphik2_curve, int verbose);

