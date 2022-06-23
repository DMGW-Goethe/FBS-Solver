import pyfbs

eos = pyfbs.PyEoStable("/media/data/Documents/PhD/B09/B09/EOS_tables/eos_HS_DD2_with_electrons.beta")

mu=1.
lamb = 0.

myFBS = pyfbs.PyFermionBosonStar(eos, 1., lamb)

rho_c = 0.0005
phi_c = 1e-6
myFBS.set_initial_conditions(rho_c, phi_c)
myFBS.bisection(1., 10., max_step=200)
myFBS.evaluate_model()
print(myFBS.get())

