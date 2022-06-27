from matplotlib import scale
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import matplotlib.tri as tri # to create a mask for a concarve shape

def scatter_plotter(sol_array):
	# extract wanted solution from the array:
	allRadii = sol_array[:,4]
	allMasses = sol_array[:,0]

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

# plot the stability region by interpolaring between all stable points:
# input an array with filtered stars (only the ones that are stable)
def plot_interpolate_stability_region(sol_filtered, z_index, maxR, maxM, triang_cutoff):

	# data coordinates and values
	allRadii = sol_filtered[:,4]
	allMasses = sol_filtered[:,0]
	allZvalues = sol_filtered[:,z_index]

	# target grid to interpolate to
	xi = np.arange(0.0 , maxR, 0.01) # range for R
	yi = np.arange(0.0 , maxM, 0.01)	# range for M
	xi,yi = np.meshgrid(xi,yi)

	
	# triangulate the data
	mytriang = tri.Triangulation(allRadii, allMasses)

	# set mask if we want to exclde some region (because our region is convex!))
	# apply mask to exclude the triangles outside of the wanted region (by essentially setting a triangle cutoff)
	apply_mask(mytriang, allRadii, allMasses, triang_cutoff)

	# plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	# interpolate
	plt.tricontourf(mytriang, allZvalues, cmap='rainbow')

	plt.xlim([0.0, maxR])
	plt.ylim([0.0, maxM])

	plt.plot(allRadii, allMasses,'k.')
	plt.title("stability region DD2", fontsize=16)
	plt.xlabel(r'$R$',fontsize=16)
	plt.ylabel(r'$M$',fontsize=16)
	cb1 = plt.colorbar(orientation="vertical")
	#cb1.ax.tick_params(labelsize=10, fontsize=16)
	cb1.set_label(label=r"$\phi_c$",size=16)#, weight='bold')

	plt.show()