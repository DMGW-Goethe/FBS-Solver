import load_data
import plotting_functions as pfuncts
import stability_curve_calc as scc

if __name__ == "__main__":

	NstarsRho = 40
	NstarsPhi = 40

	filename1 = "../plots/NbNf_test1.txt"
	filename2 = "../plots/DD2_stab_curve_test1.txt"
	filename3 = "../plots/polytrope_stab_curve_test5.txt"
	df = load_data.load_MRPhi_data(filename3)
	#print(df)
	#pfuncts.scatter_plotter(df)

	#pfuncts.plot_interpolate_stability_region(df, 2, 20., 2.5, 0.45)

	stab_curve = scc.calc_stability_curve(df, NstarsRho, NstarsPhi, 2)

	pfuncts.plot_rho_phi_stability_curve_diagram(df, stab_curve, NstarsRho, NstarsPhi)

	pfuncts.plot_interpolate_stability_region(df, 2, 20.0 , 2.5 , 0.55)

	exit()