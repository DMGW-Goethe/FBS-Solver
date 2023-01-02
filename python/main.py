import load_data
import plotting_functions as pfuncts
import stability_curve_calc as scc
import numpy as np


import matplotlib.pyplot as plt  # for testing purposes. remove later

if __name__ == "__main__":

#	filename1 = "../plots/NbNf_test1.txt"
	filename2 = "DD2_200x200_mu1_lam0.txt"
#	filename3 = "../plots/polytrope_stab_curve_test5.txt"
	df, indices = load_data.load_MRPhi_data(filename2)
	#print(df)
	#pfuncts.scatter_plotter(df)
	NstarsRho = len(np.unique(df[:,indices['rho_0']]))
	NstarsPhi = len(np.unique(df[:,indices['phi_0']]))

	#pfuncts.plot_interpolate_stability_region(df, 2, 20., 2.5, 0.45)

	print("calcing stab curve")
	stab_curve = scc.calc_stability_curve(df, indices, NstarsRho, NstarsPhi, 2)

	# colour and corresponding levels and contour lines for 2D Mass-plot:
	contourlines_colour = ['purple', 'brown', 'red', 'green', 'orange']
	contourlines_levels = [0.62, 1.2, 1.6, 2.0, 2.3]
	my_plot_filename = "test_rhoc_phic_plot.png"
	plot_title = "DD2 EOS blabla"

	rhoplotmin = df[:,indices['rho_0']][0]
	rhoplotmax = df[:,indices['rho_0']][-1]
	phiplotmin = df[:,indices['phi_0']][0]
	phiplotmax = df[:,indices['phi_0']][-1]

	print("making first plot")
	pfuncts.plot_rho_phi_stability_curve_diagram(df, indices, stab_curve, NstarsRho, NstarsPhi, contourlines_colour, contourlines_levels, rhoplotmin, phiplotmin, rhoplotmax, phiplotmax, plot_title, my_plot_filename)

	
	my_plot_filename = "test_MR_plot_filteredtest.png"
	plot_title = "stability region DD2 filtered test"
	print("making second plot")
	pfuncts.plot_interpolate_stability_region(df, indices, 2, 20.0, 2.5, 0.55, False, plot_title, my_plot_filename)

	print("filtering data")
	filtered_data_arr = scc.filter_stab_curve_data(df, indices, stab_curve)

	#print(type(df))
	#print(type(filtered_data_arr))
	print("making last plot")
	pfuncts.plot_interpolate_stability_region(filtered_data_arr, indices, 2, 20.0, 2.5, 0.55, True, plot_title, my_plot_filename)

	print("idk man")
	my_stabcurveMR = scc.stab_curve_to_MR(df, indices, stab_curve, NstarsRho, NstarsPhi)

	#print(my_stabcurveMR)

	##----------
	# (for testing purposes. remove later):
	plt.plot(my_stabcurveMR[:,1], my_stabcurveMR[:,0], c='green', linewidth=2)
	plt.show()
	##----------

	exit()
