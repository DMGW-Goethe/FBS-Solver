import load_data
import plotting_functions as pfuncts
import stability_curve_calc as scc


import matplotlib.pyplot as plt

if __name__ == "__main__":


	filename1 = "../plots/NbNf_test1.txt"
	filename2 = "../plots/DD2_MR_MRphi-plot1.txt"
	filename3 = "../plots/polytrope_stab_curve_test5.txt"
	df = load_data.load_MRPhi_data(filename3)
	#print(df)
	#pfuncts.scatter_plotter(df)

	#pfuncts.plot_interpolate_stability_region(df, 2, 20., 2.5, 0.45)

	stab_curve = scc.calc_stability_curve(df, 40, 40, 2)

	#print(stab_curve)
	Mdata = []
	Rdata = []
	Rhocdata = []
	Phicdata = []
	for i in range(len(stab_curve)):
		Mdata.append(stab_curve[i][1])
		Rdata.append(stab_curve[i][0])
		Rhocdata.append(stab_curve[i][2])
		Phicdata.append(stab_curve[i][3])

	print(Mdata)
	#plt.xlim([0,20])
	#plt.ylim([0,2.5])
	#plt.scatter(Rdata, Mdata)
	#plt.scatter(Rhocdata,Phicdata)
	#plt.show()

	exit()