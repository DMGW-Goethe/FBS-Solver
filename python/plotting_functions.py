#from matplotlib import scale
import matplotlib.pyplot as plt
import numpy as np
#from scipy.interpolate import griddata
import matplotlib.tri as tri # to create a mask for a concarve shape

import stability_curve_calc as scc # import custom python file which computes the stability curve

# plot a number of points in a scatter plot
# !unfinished!
def scatter_plotter(sol_array, indices):
	# extract wanted solution from the array:
	allRadii = sol_array[:,indices['R_F_0']]
	allMasses = sol_array[:,indices['M_T']]

	#print("\n")
	#print(allRadii)
	#print(allMasses)
	plt.scatter(allRadii, allMasses, label="FBS Configurations")
	plt.show()


# helping function to plot a concarve surface (mask out some regions)
def apply_mask(triang, x, y, triang_cutoff=0.5):
    # Mask triangles with sidelength bigger some alpha
    triangles = triang.triangles
    # Mask off unwanted triangles.
    xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
    ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
    maxi = np.max(np.sqrt((0.2*xtri)**2 + ytri**2), axis=1)
    # apply masking
    triang.set_mask(maxi > triang_cutoff)

# plot the stability region in an MR-diagram (!) by interpolaring between all stable points:
# input an array with filtered stars (only the ones that are stable)
def plot_interpolate_stability_region(sol_filtered, indices, z_index, maxR, maxM, triang_cutoff, BOOLplotpoints, myfigtitle, myfilename):

	# data coordinates and values
	allRadii = sol_filtered[:,indices['R_F_0']]
	allMasses = sol_filtered[:,indices['M_T']]
	allZvalues = sol_filtered[:,z_index]

	# target grid to interpolate to
	xi = np.arange(0.0 , maxR, 0.01) # range for R
	yi = np.arange(0.0 , maxM, 0.01) # range for M
	xi,yi = np.meshgrid(xi,yi)

	
	# triangulate the data
	mytriang = tri.Triangulation(allRadii, allMasses)

	# set mask if we want to exclde some region (because our region is convex!))
	# apply mask to exclude the triangles outside of the wanted region (by essentially setting a triangle cutoff)
	apply_mask(mytriang, allRadii, allMasses, triang_cutoff)

	# plot
	fig = plt.figure(figsize=(8,6))
	ax = fig.add_subplot(111)
	# interpolate
	plt.tricontourf(mytriang, allZvalues, cmap='rainbow')

	plt.xlim([0.0, maxR])
	plt.ylim([0.0, maxM])

	if BOOLplotpoints:
		plt.plot(allRadii, allMasses,'k.')

	plt.title(myfigtitle, fontsize=16)
	plt.xlabel(r'$R$ [km]',fontsize=16)
	plt.ylabel(r"Total Gravitational Mass [M$_\odot$]",fontsize=16)
	cb1 = plt.colorbar(orientation="vertical")
	#cb1.ax.tick_params(labelsize=10, fontsize=16)
	cb1.set_label(label=r"$\phi_c$ [Code Units]",size=16)#, weight='bold')

	#plt.show()
	plt.savefig(fname=myfilename, dpi=300)


# plot the stability curve and contour lines of constanf M in a rho_c phi_c diagram:
def plot_rho_phi_stability_curve_diagram(sol_array, indices, mystability_curve, num_rho_stars, num_phi_stars, cont_colour, cont_levels, rhoplotmax, phiplotmax, myplottitle, myfilename):

	# ------------------------------------------------------------
	# data pre-processing:

	# re-arrange the solution array:
	# create empty 3D solution array:
	ordered_sol = scc.create_3D_array(num_rho_stars, num_phi_stars, len(sol_array[0]))
	# fill the "ordered" (in a sense that it is nw a 2D array and not a simple list) solution array:
	for i in range(num_rho_stars):
		for j in range(num_phi_stars):
			for k in range(len(sol_array[0]) ):
				tmp = sol_array[i+ num_rho_stars*j][k]
				ordered_sol[i][j][k] = tmp


	M_array = scc.create_3D_array(num_rho_stars, num_phi_stars, 1) # is actually a 2D array
	for i in range(num_rho_stars):
		for j in range(num_phi_stars):
			M_array[j][i] = ordered_sol[i][j][indices['M_T']]

	# init the grid on which the 2D plot is oriented on
	rhogrid = []
	phigrid = []
	# create a grid in rho_c and phi_c using the available data:
	for i in range(num_rho_stars):
		rhogrid.append(ordered_sol[i][0][indices['rho_0']])
	for j in range(num_phi_stars):
		phigrid.append(ordered_sol[0][j][indices['phi_0']])

	# ------------------------------------------------------------
	# start with the plotting:
	plt.figure(figsize=(8,6))

	# plotting contour lines of constant M:
	contours_M = plt.contour(rhogrid, phigrid, M_array, colors=cont_colour, levels=cont_levels)

	# plot a 2D plot:
	plt.imshow(M_array, extent=[0.0, rhoplotmax, 0.0, phiplotmax], origin='lower', cmap='turbo', aspect='auto', interpolation='none', alpha=0.8)

	# plot the stability curve:
	plt.plot(mystability_curve[:,0], mystability_curve[:,1], c='black', linewidth=2)

	# add labels and stuff:
	plt.clabel(contours_M, inline=True, fontsize = 10)
	plt.xlabel(r"$\rho_c$ [Code Units]", fontsize = 15)
	plt.ylabel(r"$\phi_c$ [Code Units]", fontsize = 15)

	plt.title(myplottitle, fontsize = 15)

	# add color bar on the side of the plot
	cbar = plt.colorbar()
	cbar.set_label(r"Total Gravitational Mass [M$_\odot$]", rotation=270, fontsize=15)
	cbar.ax.get_yaxis().labelpad = 15

	#plt.show()
	plt.savefig(fname=myfilename, dpi=300)
