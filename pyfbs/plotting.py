#from matplotlib import scale
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#from scipy.interpolate import griddata
import matplotlib.tri as tri # to create a mask for a concarve shape
from matplotlib import lines
from matplotlib import cm
from matplotlib.offsetbox import AnchoredText
import pandas as pd

from . import stability_curve as scc # import custom python file which computes the stability curve
from . import data
from . import errorplot as errp

# Function that returns all solutions in df that have condition(df, indices) == True
def searchData(df, indices, condition):
	result = []

	for sol in df:
		if(sol[indices["phi_0"]] == 0 and sol[indices["rho_0"]] == 0):
			continue

		if(condition(sol, indices)):
			result.append(sol)

	return np.array(result)

# Using the above function to find all solutions that correspond to pure neutron stars
def searchPureNS(df, indices):
	def condition(sol, indices):
		if(sol[indices["phi_0"]] == 0):
			return True
		return False

	return searchData(df, indices, condition)

# Using the searchData function to find all solutions that correspond to pure boson stars
def searchPureBS(df, indices):
	def condition(sol, indices):
		if(sol[indices["rho_0"]] == 0):
			return True
		return False
		
	return searchData(df, indices, condition)


# Function that plots multiple dfs in a grid. Parameters are:
# scatterFunc -> specifies what has to be plotted.
# figHeight and figWidth specify the overall height and width. The default values seem to work nice for quadratic grids
# nrows and ncols specify the amount of rows and cols in the grid
# s -> specifies the circle size of the dots in the scatter plot
# plotPureNS -> If yes, then the line corresponding to pure neutron star solutions is also shown
# plotPureBS -> If yes, then the line corresponding to pure boson star solutions is also shown
# filterData -> If yes, then the provided data is filtered according the stability curve which is found by using scc.calc_stablity_curve
# cmap -> cmap that is used for the scatter plot
# xlim and ylim define the axis limits for all individual cells
# tickFontSize defines the font size of all ticks and the fontsize of the inset labels
def grid_Scatter(dfs, indices, scatterFunc, ylabel="", xlabel="", figHeight = None, figWidth = None, nrows = 2, ncols = 2, stabCurves=None,
                       overlay_info=None, overlay_loc=1, addColorbar = True, lockAxes = True, addLegend = True, tickFontSize=15, cmap='viridis',
                       **kwargs):
	#s = 0.3, plotPureNS = True, plotPureBS = False, filterData = True, cmap = "viridis", xlim = [0, 2.7], ylim = [1, 5e6], tickFontSize = 13, 
	# create a grid of plots
	fig = plt.figure()
	gs = fig.add_gridspec(nrows, ncols, hspace=0, wspace=0)
	if lockAxes:
		axs = gs.subplots(sharex='col', sharey='row')
	else:
		axs = gs.subplots()
	 
	# if no figHeight and width are specified, then use these values that seem to work nice for quadratic grids
	# if the grid is not quadratic, then try to adjust the default width and height to make it look nicer
	if figHeight == None:
		figHeight = 13
		if nrows > ncols:
		    figHeight /= ncols
		    figHeight *= nrows
	if figWidth == None:
		figWidth = figHeight * 1.3
		if nrows > ncols:
		    figWidth /= nrows
		    figWidth *= ncols


	fig.set_figheight(figHeight)
	fig.set_figwidth(figWidth)

	# iterate over all dfs and fill the individual cells of the grid according the provided plotFunc
	for i, ax in enumerate(axs.flatten()):
		scatterFunc(dfs[i], indices, **kwargs, pltAxs=ax, cmap=cmap, stabCurve=stabCurves[i] if not stabCurves is None else None)
		ax.tick_params(axis='both', which='major', labelsize=tickFontSize)
		
		if not overlay_info is None: #and i % np.shape(axs)[0] == 0:
		    text = ''
		    for k,v in zip(overlay_info.keys(), overlay_info.values()):
		        text += f"${v} = {(dfs[i][0,indices[k]]):.0f} $\n"
		    
		    text_box = AnchoredText(text.strip(), frameon=True, loc=overlay_loc, pad=0.1, prop={'fontsize':18}, alpha=1.)
		    plt.setp(text_box.patch, boxstyle='round', edgecolor='black', facecolor=(0., 0., 0., 0.))
		    ax.add_artist(text_box)
		    #ax.text(0.1, 0.1, text.strip(), transform=ax.transAxes, alpha=0.5, bbox=dict(boxstyle='round', pad=0.2, edgecolor='black', facecolor=(0., 0., 0., 0.)), fontsize=18)

	fig.add_subplot(111, frame_on=False)
	plt.tick_params(labelcolor="none", bottom=False, left=False)

	# add labels
	plt.xlabel(xlabel, fontsize = 22, labelpad=18)
	plt.ylabel(ylabel, fontsize = 22, labelpad=18)


	# add the color bar to the top of the grid
	if addColorbar:
		cb_ax = fig.add_axes([0.15, 0.91, 0.5, 0.014])
		norm = mpl.colors.Normalize(vmin=0,vmax=1)
		cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cb_ax, cmap = "jet", orientation = "horizontal")
		cbar.set_label(r"Dark Matter Mass Fraction", rotation=0, fontsize=22)
		cbar.ax.get_yaxis().labelpad = 20
		#cbar.ax.xaxis.set_ticks_position('top')
		cbar.ax.xaxis.set_label_position('top')
		cbar.set_ticks(np.linspace(0, 1.0, 6, endpoint=True))
		cbar.ax.tick_params(axis='both', which='major', labelsize=tickFontSize)
	
	if addLegend:
		hand, lab = axs[0,0].get_legend_handles_labels() if nrows > 1 else axs[0].get_legend_handles_labels()
		cb_ax.legend(hand, lab, loc=(1.05, -1.5), fontsize=tickFontSize)
	
	return axs

# makes a scatter plot of two indices X and y. Parameters are:
# stabCurve -> stablity curve for the given df. If none is provided then it is automatically calculated
# filterData -> If yes, then the provided data is filtered according the stability curve which is found by using scc.calc_stablity_curve
# s -> specifies the circle size of the dots in the scatter plot
# plotPureNS -> If yes, then the line corresponding to pure neutron star solutions is also shown
# plotPureBS -> If yes, then the line corresponding to pure boson star solutions is also shown
# cmap -> cmap that is used for the scatter plot
# xlim and ylim define the axis limits for all individual cells
# pltAxs -> specifies the axs object in which the data is to be plotted. If none is provided, then plt is being used
# tickFontSize defines the font size of all ticks and the fontsize of the inset labels
def scatter_XY(df, indices, X, Y, filterData = True, stabCurve = None, s = 0.3, plotPureNS = True, plotPureBS = True, cmap = "viridis", clim = [0, 1], ylim = [1, 1e6], xlim = [0, 3], pltAxs = plt, tickFontSize = 13, xlabel="", ylabel="", yscale='linear', xscale='linear'):
	# calculate the stability curve and filter the data, if needed
	if stabCurve is None and filterData:
		stabCurve = scc.calc_stability_curve(df, indices)

	filteredData = df
	if filterData:
		filteredData = scc.filter_stab_curve_data(df, indices, stabCurve)
	
	if isinstance(cmap, str):
		cmap = cm.get_cmap(cmap)

	# extract some useful information from the data file
	allNbs = filteredData[:,indices['N_B']]
	allNfs = filteredData[:,indices['N_F']]
	allRelativMassFractions = allNbs / (allNbs + allNfs)

	allX = filteredData[:,indices[X]]
	allY = filteredData[:,indices[Y]]

	pltAxs.scatter(allX, allY, s = s, alpha=1, c=allRelativMassFractions, cmap=cmap)
	pltAxs.scatter([], [], s = 20, color = cmap(0.), label="FBS Configurations") # I have to add this empty scatter plot, since otherwise no dot is being shown in the inset label for the "FBS Configurations"
	
	# if pltAxs == plt, then finish up the plot with all labels and make it look nice
	if pltAxs == plt:
		pltAxs.clim(clim)

		cbar = pltAxs.colorbar()
		cbar.set_label(r"Dark Matter Mass Fraction", rotation=270, fontsize=tickFontSize)
		cbar.ax.get_yaxis().labelpad = 15
		cbar.ax.tick_params(axis='both', which='major', labelsize=tickFontSize)
		pltAxs.xlim(xlim)
		pltAxs.ylim(ylim)
		pltAxs.yscale(yscale)
		pltAxs.xscale(xscale)

		pltAxs.ylabel(ylabel, fontsize = tickFontSize)
		pltAxs.xlabel(xlabel, fontsize = tickFontSize)
		pltAxs.tick_params(axis='both', which='major', labelsize=tickFontSize)

	# if pltAxs != plt, then only set the limits and scale, leave the rest to the user
	else:
		pltAxs.set_xlim(xlim)
		pltAxs.set_ylim(ylim)
		pltAxs.set_yscale(yscale)
		pltAxs.set_xscale(xscale)

	pltAxs.grid(alpha=0.2, linestyle="-")
	#plt.text(10, 1.5, 'preliminary', fontsize=40, color='gray', alpha=0.05,	ha='center', va='center', rotation=30)

	# plot pure quantities (pure neutron star and/or pure boson star) if required
	if plotPureNS:
		pureNS = searchPureNS(filteredData, indices)

		pureNSX = pureNS[:,indices[X]]
		pureNSY = pureNS[:,indices[Y]]
		pltAxs.plot(pureNSX, pureNSY, label = "Pure DD2", color = cmap(0.), linestyle = "-", linewidth = 2)
		
	if plotPureBS:
		pureBS = searchPureBS(filteredData, indices)

		pureBSX = pureBS[:,indices[X]]
		pureBSY = pureBS[:,indices[Y]]
		pltAxs.plot(pureBSX, pureBSY, label = "Pure BS", color = cmap(1.), linestyle = "-", linewidth = 2, alpha=0.5)
	
	# add the final plot label (has to happen at the end because only now we know whether pure ns and/or bs was plotted)
	if pltAxs == plt:
		pltAxs.legend(loc = "upper right", fontsize = tickFontSize)


def scatter_Tidal(df, indices, filterData = True, stabCurve = None, s = 0.3, plotPureNS = True, plotPureBS = True, cmap = "viridis", clim = [0, 1], ylim = [1, 1e6], xlim = [0, 3], pltAxs = plt, tickFontSize = 13):
	df, indices = data.add_Lambda_tidal(df, indices)
	
	scatter_XY(df, indices, 'M_T', 'Lambda_tidal', filterData=filterData, stabCurve=stabCurve, s=s, plotPureNS=plotPureNS, plotPureBS=plotPureBS, 
				cmap=cmap, clim=clim, ylim=ylim, xlim=xlim, pltAxs=pltAxs, tickFontSize=tickFontSize, xlabel=r"Total Gravitational Mass [M$_\odot$]", ylabel=r"Tidal Deformability $\Lambda$", yscale='log')
				

def scatter_MR(df, indices, filterData = True, stabCurve = None, s = 0.3, plotPureNS = True, cmap = "viridis", clim = [0, 1], ylim = [1, 3], xlim = [0, 30], pltAxs = plt, tickFontSize = 13, plotPureBS=False):
	scatter_XY(df, indices, 'R_F', 'M_T', filterData=filterData, stabCurve=stabCurve, s=s, plotPureNS=plotPureNS, plotPureBS=False, 
				cmap=cmap, clim=clim, ylim=ylim, xlim=xlim, pltAxs=pltAxs, tickFontSize=tickFontSize, xlabel=r"Fermionic Radius [km]", ylabel=r"Total Gravitational Mass [M$_\odot$]")


def scatter_MRgrav(df, indices, filterData = True, stabCurve = None, s = 0.3, plotPureNS = True, plotPureBS=True, cmap = "viridis", clim = [0, 1], ylim = [1, 3], xlim = [0, 30], pltAxs = plt, tickFontSize = 13):
	df, indices = data.add_R_max(df, indices)
	
	scatter_XY(df, indices, 'R_m', 'M_T', filterData=filterData, stabCurve=stabCurve, s=s, plotPureNS=plotPureNS, plotPureBS=plotPureBS, 
				cmap=cmap, clim=clim, ylim=ylim, xlim=xlim, pltAxs=pltAxs, tickFontSize=tickFontSize, xlabel=r"Gravitational Radius [km]", ylabel=r"Total Gravitational Mass [M$_\odot$]")


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

	# set mask if we want to exclde some region (because our region is not convex!))
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


# plots a density plot in the rho-phi plane of the quantity specified with index
# df and indices are the data file and the indices that relate the entries to physical data
# index is the index which is to be plotted in the density plot
# cmap is the color map
# clim specifies the min and max values to  be colored
# scale specifies whether linear or log scaling is to be used
# stabCurve is the stabilty curve. It gets calculated automatically if needed
# plotStabCurve -> whether or not the stability curve is to be plotted
# cBarLabel -> label for the color bar
# autoContour -> if yes, then the levels for the contour lines are chosen automatically
# contourLevels -> manually specify the contour line levels
# contourColors -> manually specify the contour line colors
# tickFontSize specifies the size of the tick labels
# xlim and ylim specify the rho and phi range, respetively
def plotRhoPhi(df, indices, index, cmap = "viridis", clim = None, scale = "linear", stabCurve = None, plotStabCurve = True, cbarLabel = "", autoContour = False, contourColors = ['green', 'brown', '#EE2C2C', "#00688B", "#E8E8E8"], contourLevels = None, tickFontSize = 13, xlim = None, ylim = None):
	
	# first, get some useful quantities from the list of all values
	numStarsRho = len(np.unique(df[:,indices['rho_0']]))

	rhoVals = np.unique(df[:,indices["rho_0"]])
	phiVals = np.unique(df[:,indices["phi_0"]])

	if stabCurve is None and plotStabCurve:
		stabCurve = scc.calc_stability_curve(df, indices)

	data_frame = pd.DataFrame(df[:, [indices[index], indices['rho_0'], indices['phi_0']]], columns=[index, 'rho_0', 'phi_0'])
	plotArray = data_frame.pivot_table(values=index, columns='rho_0', index='phi_0')
	
	if scale == "log":
		plotArray = np.log(plotArray)

	if contourLevels is None and autoContour:
		minVal = min(allPlotArray)
		maxVal = max(allPlotArray)
		contourLevels = np.linspace(minVal, maxVal, 7)[1:-1]

	contourPlot = plt.contour(rhoVals, phiVals, plotArray, colors=contourColors, levels=contourLevels)
	plt.clabel(contourPlot, inline=True, fontsize = 10)

	xx, yy = np.meshgrid(rhoVals, phiVals)
	plt.pcolormesh(xx, yy, plotArray, cmap = cmap)
	if not clim is None:
		plt.clim(clim)

	if not xlim is None:
		plt.xlim(xlim)
	if not ylim is None:
		plt.ylim(ylim)

	plt.xlabel(r"$\rho_c$ [Code Units]", fontsize = 15)
	plt.ylabel(r"$\phi_c$ [Code Units]", fontsize = 15)
	
	cbar = plt.colorbar()
	cbar.set_label(cbarLabel, rotation=270, fontsize=15)
	cbar.ax.get_yaxis().labelpad = 15

	if plotStabCurve:
		plt.plot(stabCurve[:,0], stabCurve[:,1], color = "black", linewidth = 2)

# plots the relative numerical errors, of two quantities, between the full system and the effective EoS system with respect to lambda_int
# filenames is a list containing all filenames of fullsys and effsys, with every other filename corresponding to the fullsys and effsys respectively
# Lambda_array contains all Lambda_int values in ascending order matching eqch pair of effsys/fullsys data points
def plot_error_comparison_effsys_fullsys(filenames, Lambda_array, index1, index2, figname = "plot_error_comparison_effsys_fullsys.png", debug = False):

	# reading in the data using the filenames and assigning indices
	number_of_files = len(filenames)
	df = [None]*number_of_files
	indices = [None]*number_of_files

	for i in range(len(filenames)):
		df[i], indices[i] = data.load_file(filenames[i])

	indices_fullsys = indices[0]
	indices_effsys = indices[1]
	# computation of the relative errors
	min_err1, max_err1, median_err1, one_sigma_err1, two_sigma_err1 = errp.calc_error_curves(df, indices_fullsys, indices_effsys, index1)
	min_err2, max_err2, median_err2, one_sigma_err2, two_sigma_err2 = errp.calc_error_curves(df, indices_fullsys, indices_effsys, index2)

	# start the plotting:
	if (debug):
		print("for index1 = ", index1, ":\n max_err: \t", max_err1, ",\n min_err: \t", min_err1, ",\n median_err: \t", median_err1, ",\n 1_sigma_err: \t", one_sigma_err1, ",\n 2_sigma_err: \t", two_sigma_err1, ",\n Lambda_int: \t", Lambda_array)
		print("for index2 = ", index2, ":\n max_err: \t", max_err2, ",\n min_err: \t", min_err2, ",\n median_err: \t", median_err2, ",\n 1_sigma_err: \t", one_sigma_err2, ",\n 2_sigma_err: \t", two_sigma_err2, ",\n Lambda_int: \t", Lambda_array)

	# hardcoded parameters (might change later)
	nrows = 2
	ncols = 1
	fig = plt.figure()
	gs = fig.add_gridspec(nrows, ncols, hspace=0.1, wspace=0)
	axs = gs.subplots(sharex='col')#, sharey='row')

	# populate the first subplot
	axs[0].fill_between(Lambda_array, min_err1, max_err1, alpha = 0.15, color="blue", label = "$M_{tot}$ error region")
	axs[0].fill_between(Lambda_array, min_err1, two_sigma_err1, alpha = 0.15, color="blue")#, label = "$M_{tot}$ $2\sigma$ error region")#, hatch='x')
	axs[0].fill_between(Lambda_array, min_err1, one_sigma_err1, alpha = 0.15, color="blue") #, label = "$1\sigma$ region"
	axs[0].plot(Lambda_array, median_err1, color="blue", linewidth=1.75, linestyle ="--")#, label="median of $M_{tot}$-error")

	# populate the second subplot
	axs[1].fill_between(Lambda_array, min_err2, max_err2, alpha = 0.15, color="red", label = "$\Lambda \:\:\:\:\:\;$ error region")
	axs[1].fill_between(Lambda_array, min_err2, two_sigma_err2, alpha = 0.15, color="red")#, label = "$\Lambda \:\:\:\:\:\;$ $2\sigma$ error region")
	axs[1].fill_between(Lambda_array, min_err2, one_sigma_err2, alpha = 0.15, color="red") #, label = "$1\sigma$ region"
	axs[1].plot(Lambda_array, median_err2, color="red", linewidth=1.75, linestyle ="--")#, label="median of $\Lambda$-error")

	# set scale, plotlim, legend ...
	axs[0].set_yscale("log") # "log" or "linear"
	axs[1].set_yscale("log") # "log" or "linear"
	axs[0].set_ylim(1e-3,2e0)
	axs[1].set_ylim(1e-3,1e2)
	axs[0].set_xlim(0,300)
	axs[1].set_xlim(0,300)

	axs[0].grid(alpha=0.15, linestyle="--", c="black")
	axs[1].grid(alpha=0.15, linestyle="--", c="black")

	axs[0].legend(loc="lower left")
	axs[1].legend(loc="lower left")

	# shared x-axis label of the plot
	plt.xlabel("$\Lambda_{\mathrm{int}} = \lambda / (8 \pi m^2)$", fontsize = 14)
	# shared y-axis label
	fig.text(0.02, 0.5, "$\epsilon_{\mathrm{rel}} := ||Q_{\mathrm{full}}-Q_{\mathrm{eff}}|| / Q_{\mathrm{full}}$ , $Q = \{ M_{\mathrm{tot}}, \Lambda \}$", fontsize = 14, va='center', rotation='vertical')

	# prevent labels on the "inside" between the subplots:
	for ax in axs:
		ax.label_outer()

	# only saving pic is possible, plt.show() is not supported by matplotlib in this case
	if (debug):
		print("saving figure as: " + figname)
	plt.savefig(figname, dpi=400, bbox_inches='tight')	