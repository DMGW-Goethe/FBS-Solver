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

from matplotlib import rcParams
from matplotlib import rcParamsDefault
from matplotlib import rc

import pyfbs

rcParams.update(rcParamsDefault)

#rc('font', **{'family': 'CMU Sans Serif', 'CMU Sans Serif': ['CMUSansSerif']})
rc('font', **{'family': 'CMU Serif'})
rcParams["mathtext.fontset"] = "cm"

def latex_float(f):
    float_str = "{0:.3g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

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
def searchPureNS(df, indices, indexOffset = 0):
	allPhiVals = np.unique(df[:,indices["phi_0"]])
	desiredPhi = np.partition(allPhiVals, indexOffset)[indexOffset]

	def condition(sol, indices):
		if(sol[indices["phi_0"]] == desiredPhi):
			return True
		return False

	return searchData(df, indices, condition)

# Using the searchData function to find all solutions that correspond to pure boson stars
def searchPureBS(df, indices, indexOffset = 0):
	allRhoVals = np.unique(df[:,indices["rho_0"]])
	desiredRho = np.partition(allRhoVals, indexOffset)[indexOffset]

	def condition(sol, indices):
		if(sol[indices["rho_0"]] == desiredRho):
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
					   overlay_info=None, overlay_loc=1, addColorbar = True, lockAxes = True, addLegend = True, tickFontSize=15, cmap='viridis', clim = [0, 1], cmapRange=[0, 1],
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
		plotHandles, plotLabels = scatterFunc(dfs[i], indices, clim = clim, pltAxs=ax, cmap=cmap, cmapRange=cmapRange, stabCurve=stabCurves[i] if not stabCurves is None else None, **kwargs)
		ax.tick_params(axis='both', which='major', labelsize=tickFontSize)
		
		if not overlay_info is None: #and i % np.shape(axs)[0] == 0:
			text = ''
			for k,v in zip(overlay_info.keys(), overlay_info.values()):
				if(k is "mu"):
					text += f"${v} = {latex_float(dfs[i][0,indices[k]] * 1.34e-10)}$ eV\n"
					#text += f"${v} = {(dfs[i][0,indices[k]] * 1.34e-10):.2e} $\n"
				else:
					text += f"${v} = {(dfs[i][0,indices[k]]):.0f} $\n"
			
			text_box = AnchoredText(text.strip(), frameon=True, loc=overlay_loc, pad=0.1, prop={'fontsize':12}, alpha=1., borderpad=0.8)
			plt.setp(text_box.patch, boxstyle='round', edgecolor='black', facecolor=(1., 1., 1., 1.))
			ax.add_artist(text_box)
			#ax.text(0.1, 0.1, text.strip(), transform=ax.transAxes, alpha=0.5, bbox=dict(boxstyle='round', pad=0.2, edgecolor='black', facecolor=(0., 0., 0., 0.)), fontsize=18)

	fig.add_subplot(111, frame_on=False)
	plt.tick_params(labelcolor="none", bottom=False, left=False)

	# add labels
	plt.xlabel(xlabel, fontsize = 22, labelpad=18)
	plt.ylabel(ylabel, fontsize = 22, labelpad=18)


	# add the color bar to the top of the grid
	if addColorbar:
		cb_ax = fig.add_axes([0.15, 0.925, 0.5, 0.014])
		norm = mpl.colors.Normalize(vmin=0,vmax=clim[1])
		cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cb_ax, cmap = cmap, orientation = "horizontal")

		cbar.ax.set_xlim(0, 1)
		cbar.set_label(r"Dark Matter Mass Fraction", rotation=0, fontsize=22)
		cbar.ax.get_xaxis().labelpad = 6
		#cbar.ax.xaxis.set_ticks_position('top')
		cbar.ax.xaxis.set_label_position('top')
		cbar.set_ticks(np.linspace(0, 1.0, 6, endpoint=True))
		cbar.ax.tick_params(axis='both', which='major', labelsize=tickFontSize)
	
	if addLegend:
		#hand, lab = axs[0, 0].get_legend_handles_labels() if nrows > 1 else axs[0].get_legend_handles_labels()
		cb_ax.legend(handles = plotHandles, labels = plotLabels, loc=(1.05, -1.5), fontsize=tickFontSize)
	
	return fig, axs

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
def scatter_XY(df, indices, X, Y, rasterized=False, filterData = True, stabCurve = None, s = 0.3, plotPureNS = True, plotPureBS = True, cmap = "viridis", clim = [0, 1], ylim = [1, 1e6], xlim = [0, 3], pltAxs = plt, tickFontSize = 13, xlabel="", ylabel="", yscale='linear', xscale='linear', cmapRange=[0, 1], pureNSIndexOffset = 0, legendloc = "upper left", showBuchdahlLimit = False):
	# calculate the stability curve and filter the data, if needed
	if stabCurve is None and filterData:
		stabCurve = scc.calc_stability_curve(df, indices)

	filteredData = df
	if filterData:
		filteredData = scc.filter_stab_curve_data(df, indices, stabCurve)
	
	if isinstance(cmap, str):
		cmap = cm.get_cmap(cmap)

	plotHandles = []
	plotLables = []

	# extract some useful information from the data file
	allNbs = filteredData[:,indices['N_B']]
	allNfs = filteredData[:,indices['N_F']]
	allRelativMassFractions = allNbs / (allNbs + allNfs)

	allX = filteredData[:,indices[X]]
	allY = filteredData[:,indices[Y]]

	pltAxs.scatter(allX, allY, s = s, alpha=1, c=allRelativMassFractions, cmap=cmap, vmin=clim[0], vmax=clim[1], rasterized=rasterized)

	# if pltAxs == plt, then finish up the plot with all labels and make it look nice
	if pltAxs == plt:
		pltAxs.clim(clim)

		cbar = pltAxs.colorbar()
		cbar.ax.set_ylim(cmapRange[0], cmapRange[1])
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

	if(showBuchdahlLimit):
		p1 = xlim
		p2 = [xlim[0]*4/9, xlim[1]*4/9]
		pltAxs.fill_between(x=p1, y1=p2, y2=ylim[1], color="black", alpha=0.1)
		#plt.fill_between([0 xlim[1]], y1=[0, 2], y2=[3,3])
		pltAxs.plot(p1, p2, color="black")

	scatterHandle = pltAxs.scatter([], [], s = 20, color = "black", label="FBS Configurations", cmap=cmap) # I have to add this empty scatter plot, since otherwise no dot is being shown in the inset label for the "FBS Configurations"
	plotHandles.append(scatterHandle)
	plotLables.append("FBS Configurations")

	pltAxs.grid(alpha=0.2, linestyle="--")
	#plt.text(10, 1.5, 'preliminary', fontsize=40, color='gray', alpha=0.05,	ha='center', va='center', rotation=30)

	# plot pure quantities (pure neutron star and/or pure boson star) if required
	if plotPureBS:
		pureBS = searchPureBS(filteredData, indices)

		pureBS = np.array([sol for sol in pureBS if np.abs(sol[indices[Y]]) < 10 * ylim[1] and np.abs(sol[indices[Y]]) > ylim[0]])
		pureBSX = pureBS[:,indices[X]]
		pureBSY = pureBS[:,indices[Y]]

		pltAxs.plot(pureBSX, pureBSY, label = "Pure BS", color = cmap(cmapRange[1]/clim[1]), linestyle = "-", linewidth = 4, alpha=1)
		pltAxs.plot(pureBSX, pureBSY, color = "white", linestyle = "-", linewidth = 1.5)

		pureBSlineLegend1 = lines.Line2D([], [], linewidth=4, linestyle="-", color=cmap(cmapRange[1]/clim[1]))
		pureBSlineLegend2 = lines.Line2D([], [], linewidth=1.5, linestyle="-", color='white')

		plotHandles.append((pureBSlineLegend1, pureBSlineLegend2))
		plotLables.append("Pure BS")

	if plotPureNS:
		pureNS = searchPureNS(filteredData, indices, pureNSIndexOffset)
		pureNS = np.array([sol for sol in pureNS if np.abs(sol[indices[Y]]) < 10 * ylim[1] and np.abs(sol[indices[Y]]) > ylim[0] and np.abs(sol[indices[X]]) > 1e-9 and np.abs(sol[indices[X]]) < 10 * xlim[1]])
		
		pureNSX = pureNS[:,indices[X]]
		pureNSY = pureNS[:,indices[Y]]

		#print(pureNSX)
		pureDD2handle, = pltAxs.plot(pureNSX, pureNSY, label = "Pure DD2", color = cmap(0.), linestyle = "-", linewidth = 4)#, gapcolor="yellow", dashes=[3, 3])
		pltAxs.plot(pureNSX, pureNSY, color = "white", linestyle = "-", linewidth = 1.5)
		
		pureNSlineLegend1 = lines.Line2D([], [], linewidth=4, linestyle="-", color=cmap(0))
		pureNSlineLegend2 = lines.Line2D([], [], linewidth=1.5, linestyle="-", color='white')

		plotHandles.append((pureNSlineLegend1, pureNSlineLegend2))
		plotLables.append("Pure DD2")

	# add the final plot label (has to happen at the end because only now we know whether pure ns and/or bs was plotted)
	if pltAxs == plt:
		pltAxs.legend(handles = plotHandles, labels = plotLables, loc = legendloc, fontsize = 12, framealpha=1)
	
	return plotHandles, plotLables


def scatter_Tidal(df, indices, **kwargs):
	df, indices = data.add_Lambda_tidal(df, indices)
	
	return scatter_XY(df, indices, 'M_T', 'Lambda_tidal', xlabel=r"Total Gravitational Mass [M$_\odot$]", ylabel=r"Tidal Deformability $\Lambda$", yscale='log', **kwargs)
				
'''
def scatter_MR(df, indices, filterData = True, stabCurve = None, s = 0.3, plotPureNS = True, cmap = "viridis", clim = [0, 1], ylim = [1, 3], xlim = [0, 30], pltAxs = plt, tickFontSize = 13, plotPureBS=False, **kwargs):
	scatter_XY(df, indices, 'R_F', 'M_T', filterData=filterData, stabCurve=stabCurve, s=s, plotPureNS=plotPureNS, plotPureBS=False, 
				cmap=cmap, clim=clim, ylim=ylim, xlim=xlim, pltAxs=pltAxs, tickFontSize=tickFontSize, xlabel=r"Fermionic Radius [km]", ylabel=r"Total Gravitational Mass [M$_\odot$]", **kwargs)
'''

def scatter_MR(df, indices, **kwargs):
	return scatter_XY(df, indices, 'R_F', 'M_T', xlabel=r"Fermionic Radius [km]", ylabel=r"Total Gravitational Mass [M$_\odot$]", **kwargs)


def scatter_MRgrav(df, indices, **kwargs):
	df, indices = data.add_R_max(df, indices)
	
	return scatter_XY(df, indices, 'R_m', 'M_T', xlabel=r"Gravitational Radius [km]", ylabel=r"Total Gravitational Mass [M$_\odot$]", **kwargs)


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
def plotRhoPhi(df, indices, index, xlabel = r"$\rho_c$ [Code Units]", ylabel = r"$\phi_c$ [Code Units]", cmap = "plasma", clim = None, scale = "linear", stabCurve = None, plotStabCurve = True, cbarLabel = "", autoContour = False, contourColors = ['green', 'brown', '#EE2C2C', "#00688B", "#E8E8E8"], contourLevels = None, tickFontSize = 13, xlim = None, ylim = None, manualContourPos = None):
	
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
		minVal = min(plotArray)
		maxVal = max(plotArray)
		contourLevels = np.linspace(minVal, maxVal, 7)[1:-1]

	contourPlot = plt.contour(rhoVals, phiVals, plotArray, colors=contourColors, levels=contourLevels)
	plt.clabel(contourPlot, inline=True, fontsize = 11, manual=manualContourPos)
	plt.tick_params(axis='both', which='major', labelsize=14)

	xx, yy = np.meshgrid(rhoVals, phiVals)
	plt.pcolormesh(xx, yy, plotArray, cmap = cmap, alpha = 0.8, linewidth=0, rasterized=True)
	if not clim is None:
		plt.clim(clim)

	if not xlim is None:
		plt.xlim(xlim)
	if not ylim is None:
		plt.ylim(ylim)

	plt.xlabel(xlabel, fontsize = tickFontSize)
	plt.ylabel(ylabel, fontsize = tickFontSize)
	
	cbar = plt.colorbar()
	cbar.set_label(cbarLabel, rotation=270, fontsize=tickFontSize)
	cbar.ax.get_yaxis().labelpad = 20
	cbar.ax.tick_params(axis='both', which='major', labelsize=14)

	if plotStabCurve:
		plt.plot(stabCurve[:,0], stabCurve[:,1], color = "black", linewidth = 2)
