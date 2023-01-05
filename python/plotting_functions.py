#from matplotlib import scale
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#from scipy.interpolate import griddata
import matplotlib.tri as tri # to create a mask for a concarve shape
from matplotlib import lines

import stability_curve_calc as scc # import custom python file which computes the stability curve


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
# plotFunc -> specifies what has to be plotted. Currently pfuncts.plotTidal and pfuncts.plotMR are supported
# figHeight and figWidth specify the overall height and width. The default values seem to work nice for quadratic grids
# nrows and ncols specify the amount of rows and cols in the grid
# s -> specifies the circle size of the dots in the scatter plot
# plotPureNS -> If yes, then the line corresponding to pure neutron star solutions is also shown
# plotPureBS -> If yes, then the line corresponding to pure boson star solutions is also shown
# filterData -> If yes, then the provided data is filtered according the stability curve which is found by using scc.calc_stablity_curve
# cmap -> cmap that is used for the scatter plot
# xlim and ylim define the axis limits for all individual cells
# tickFontSize defines the font size of all ticks and the fontsize of the inset labels
def plotGrid(dfs, indices, plotFunc, figHeight = None, figWidth = None, nrows = 2, ncols = 2, s = 0.3, plotPureNS = True, plotPureBS = False, filterData = True, cmap = "viridis", xlim = [0, 2.7], ylim = [1, 5e6], tickFontSize = 13):
	
	# create a grid of plots
	fig = plt.figure()
	gs = fig.add_gridspec(nrows, ncols, hspace=0, wspace=0)
	axs = gs.subplots(sharex=True, sharey=True)
	#fig.suptitle('Sharing both axes')
	 
	# if no figHeight and width are specified, then use these values that seem to work nice for quadratic grids
	if figHeight == None:
		figHeight = 13
	if figWidth == None:
		figWidth = figHeight * 1.3

	# if the grid is not quadratic, then try to adjust the default width and height to make it look nicer
	if nrows > ncols:
		figHeight /= ncols
		figHeight *= nrows

	if nrows > ncols:
		figWidth /= nrows
		figWidth *= ncols

	fig.set_figheight(figHeight)
	fig.set_figwidth(figWidth)

	# iterate over all dfs and fill the individual cells of the grid according the provided plotFunc
	for i, ax in enumerate(axs.flatten()):
		plotFunc(dfs[i], indices, filterData = filterData, pltAxs = ax, cmap=cmap, xlim=xlim, ylim=ylim, s=s, plotPureBS=plotPureBS, plotPureNS=plotPureNS)
		ax.tick_params(axis='both', which='major', labelsize=tickFontSize)

	fig.add_subplot(111, frame_on=False)
	plt.tick_params(labelcolor="none", bottom=False, left=False)

	# add labels
	plt.xlabel(r"Total Gravitational Mass [M$_\odot$]", fontsize = 22, labelpad=18)
	plt.ylabel(r"Tidal Deformability $\Lambda$", fontsize = 22, labelpad=18)

	# the inset label should be in the upper right cell (I just think that looks nice)
	if ncols == 1 and nrows != 1:
		axs[0].legend(loc = "upper right", fontsize = tickFontSize)
	elif ncols != 1 and nrows == 1:
		axs[-1].legend(loc = "upper right", fontsize = tickFontSize)
	else:
		axs[0, -1].legend(loc = "upper right", fontsize = tickFontSize)

	# add the color bar to the top of the grid
	cb_ax = fig.add_axes([0.25, 0.91, 0.5, 0.014])
	norm = mpl.colors.Normalize(vmin=0,vmax=1)
	cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cb_ax, cmap = "jet", orientation = "horizontal")
	cbar.set_label(r"Dark Matter Mass Fraction", rotation=0, fontsize=22)
	cbar.ax.get_yaxis().labelpad = 20

	#cbar.ax.xaxis.set_ticks_position('top')
	cbar.ax.xaxis.set_label_position('top')
	
	cbar.set_ticks(np.linspace(0, 1.0, 6, endpoint=True))
	cbar.ax.tick_params(axis='both', which='major', labelsize=tickFontSize)


# makes a scatter plot of the dimensionless tidal deformabiility. Parameters are:
# stabCurve -> stablity curve for the given df. If none is provided then it is automatically calculated
# filterData -> If yes, then the provided data is filtered according the stability curve which is found by using scc.calc_stablity_curve
# s -> specifies the circle size of the dots in the scatter plot
# plotPureNS -> If yes, then the line corresponding to pure neutron star solutions is also shown
# plotPureBS -> If yes, then the line corresponding to pure boson star solutions is also shown
# cmap -> cmap that is used for the scatter plot
# xlim and ylim define the axis limits for all individual cells
# pltAxs -> specifies the axs object in which the data is to be plotted. If none is provided, then plt is being used
# tickFontSize defines the font size of all ticks and the fontsize of the inset labels
def plotTidal(df, indices, filterData = True, stabCurve = None, s = 0.3, plotPureNS = True, plotPureBS = False, cmap = "viridis", clim = [0, 1], ylim = [1, 1e6], xlim = [0, 3], pltAxs = plt, tickFontSize = 13):
	# calculate the stability curve and filter the data, if needed
	if stabCurve == None and filterData:
		stabCurve = scc.calc_stability_curve(df, indices)

	filteredData = df
	if filterData:
		filteredData = scc.filter_stab_curve_data(df, indices, stabCurve)

	# extract some useful information from the data file
	allNbs = filteredData[:,indices['N_B']]
	allNfs = filteredData[:,indices['N_F']]
	allRelativMassFractions = allNbs / (allNbs + allNfs)

	allRadii = filteredData[:,indices['R_F']]
	allMasses = filteredData[:,indices['M_T']]

	allTlns = filteredData[:,indices["lambda_tidal"]]
	allDimlessTlns = allTlns / allMasses**5

	pltAxs.scatter(allMasses, allDimlessTlns, s = s, alpha=1, c=allRelativMassFractions, cmap=cmap)
	pltAxs.scatter([], [], s = 20, c = "black", label="FBS Configurations") # I have to add this empty scatter plot, since otherwise no dot is being shown in the inset label for the "FBS Configurations"
	
	# if pltAxs == plt, then finish up the plot with all labels and make it look nice
	if pltAxs == plt:
		pltAxs.clim(clim)

		cbar = pltAxs.colorbar()
		cbar.set_label(r"Dark Matter Mass Fraction", rotation=270, fontsize=tickFontSize)
		cbar.ax.get_yaxis().labelpad = 15
		cbar.ax.tick_params(axis='both', which='major', labelsize=tickFontSize)
		pltAxs.xlim(xlim)
		pltAxs.ylim(ylim)
		pltAxs.yscale("log")

		pltAxs.ylabel(r"Tidal Deformability $\Lambda$", fontsize = tickFontSize)
		pltAxs.xlabel(r"Total Gravitational Mass [M$_\odot$]", fontsize = tickFontSize)
		pltAxs.tick_params(axis='both', which='major', labelsize=tickFontSize)

	# if pltAxs != plt, then only set the limits and scale, leave the rest to the user
	else:
		pltAxs.set_xlim(xlim)
		pltAxs.set_ylim(ylim)
		pltAxs.set_yscale("log")

	pltAxs.grid(alpha=0.2, linestyle="-")
	#plt.text(10, 1.5, 'preliminary', fontsize=40, color='gray', alpha=0.05,	ha='center', va='center', rotation=30)

	# plot pure quantities (pure neutron star and/or pure boson star) if required
	if plotPureNS:
		pureNS = searchPureNS(filteredData, indices)

		pureNSMasses = pureNS[:,indices["M_T"]]
		pureNSTlns = pureNS[:,indices["lambda_tidal"]]
		pltAxs.plot(pureNSMasses, pureNSTlns / pureNSMasses**5, label = "Pure DD2", c = "black", linestyle = "-", linewidth = 2)
			
	if plotPureBS:
		pureBS = searchPureBS(filteredData, indices)

		pureBSMasses = pureBS[:,indices["M_T"]]
		pureBSTlns = pureBS[:,indices["lambda_tidal"]]
		pltAxs.plot(pureBSMasses, pureBSTlns / pureBSMasses**5, label = "Pure BS", c = "darkred", linestyle = "-", linewidth = 2)
			
	# add the final plot label (has to happen at the end because only now we know whether pure ns and/or bs was plotted)
	if pltAxs == plt:
		pltAxs.legend(loc = "upper right", fontsize = tickFontSize)

	
# makes a scatter plot of the total gravitational mass and fermionic radii. Parameters are:
# stabCurve -> stablity curve for the given df. If none is provided then it is automatically calculated
# filterData -> If yes, then the provided data is filtered according the stability curve which is found by using scc.calc_stablity_curve
# s -> specifies the circle size of the dots in the scatter plot
# plotPureNS -> If yes, then the line corresponding to pure neutron star solutions is also shown
# plotPureBS -> not relevant for this plot, as the fermionic radii are not defined in this case (basically they are 0). This argument is only added to make it work better together with the plotGrid function above
# cmap -> cmap that is used for the scatter plot
# xlim and ylim define the axis limits for all individual cells
# pltAxs -> specifies the axs object in which the data is to be plotted. If none is provided, then plt is being used
# tickFontSize defines the font size of all ticks and the fontsize of the inset labels
def plotMR(df, indices, s = 0.3, filterData = True, plotPureNS = True, plotPureBS = None, stabCurve = None, xlim = [0, 20], ylim = [0, 3], clim = [0, 1], cmap = "viridis", pltAxs = plt, tickFontSize = 13):
	if stabCurve == None and filterData:
		stabCurve = scc.calc_stability_curve(df, indices)

	filteredData = df
	if filterData:
		filteredData = scc.filter_stab_curve_data(df, indices, stabCurve)

	# extract wanted solution from the array:
	allNbs = filteredData[:,indices['N_B']]
	allNfs = filteredData[:,indices['N_F']]
	allRelativMassFractions = allNbs / (allNbs + allNfs)

	allRadii = filteredData[:,indices['R_F']]
	allMasses = filteredData[:,indices['M_T']]

	pltAxs.scatter(allRadii, allMasses, s = s, alpha=1, c=(allRelativMassFractions), cmap=cmap)
	pltAxs.scatter([], [], s = 20, c = "black", label="FBS Configurations") # I have to add this empty scatter plot, since otherwise no dot is being shown in the inset label for the "FBS Configurations"

	# if pltAxs == plt, then finish up the plot with all labels and make it look nice
	if pltAxs == plt:
		pltAxs.clim(clim)

		cbar = pltAxs.colorbar()
		cbar.set_label(r"Dark Matter Mass Fraction", rotation=270, fontsize=tickFontSize)
		cbar.ax.get_yaxis().labelpad = 15
		cbar.ax.tick_params(axis='both', which='major', labelsize=tickFontSize)
		pltAxs.xlim(xlim)
		pltAxs.ylim(ylim)

		pltAxs.ylabel(r"Total Gravitational Mass [M$_\odot$]", fontsize = tickFontSize)
		pltAxs.xlabel(r"Fermionic Radius [km]", fontsize = tickFontSize)
		pltAxs.tick_params(axis='both', which='major', labelsize=tickFontSize)

	#plt.text(10, 1.5, 'preliminary', fontsize=40, color='gray', alpha=0.05,	ha='center', va='center', rotation=30)
	else:
		pltAxs.set_xlim(xlim)
		pltAxs.set_ylim(ylim)

	#plt.text(10, 1.5, 'preliminary', fontsize=40, color='gray', alpha=0.05,	ha='center', va='center', rotation=30)
	pltAxs.grid(alpha=0.2, linestyle="-")

	if plotPureNS:
		pureNS = searchPureNS(filteredData, indices)

		pureNSMasses = pureNS[:,indices["M_T"]]
		pureNSRadii = pureNS[:,indices["R_F"]]
		pltAxs.plot(pureNSRadii, pureNSMasses, label = "Pure DD2", c = "black", linestyle = "-", linewidth = 2)
		
	# add the final plot label (has to happen at the end because only now we know whether pure ns and/or bs was plotted)
	if pltAxs == plt:
		pltAxs.legend(loc = "upper left", fontsize = tickFontSize)

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

	if stabCurve == None and plotStabCurve == True:
		stabCurve = scc.calc_stability_curve(df, indices)

	allPlotArray = df[:,indices[index]]
	if scale == "log":
		allPlotArray = np.log(allPlotArray)

	plotArray = scc.refactorList(allPlotArray, numStarsRho)
	if(df[0, indices["rho_0"]] == df[1, indices["rho_0"]]):
		plotArray = np.transpose(plotArray)
	
	if contourLevels == None and autoContour == True:
		minVal = min(allPlotArray)
		maxVal = max(allPlotArray)
		contourLevels = np.linspace(minVal, maxVal, 7)[1:-1]

	contourPlot = plt.contour(rhoVals, phiVals, plotArray, colors=contourColors, levels=contourLevels)
	plt.clabel(contourPlot, inline=True, fontsize = 10)

	xx, yy = np.meshgrid(rhoVals, phiVals)
	plt.pcolormesh(xx, yy, plotArray, cmap = cmap)
	if clim != None:
		plt.clim(clim)

	if xlim != None:
		plt.xlim(xlim)
	if ylim != None:
		plt.ylim(ylim)

	plt.xlabel(r"$\rho_c$ [Code Units]", fontsize = 15)
	plt.ylabel(r"$\phi_c$ [Code Units]", fontsize = 15)
	
	cbar = plt.colorbar()
	cbar.set_label(cbarLabel, rotation=270, fontsize=15)
	cbar.ax.get_yaxis().labelpad = 15

	if plotStabCurve:
		plt.plot(stabCurve[:,0], stabCurve[:,1], color = "black", linewidth = 2)
