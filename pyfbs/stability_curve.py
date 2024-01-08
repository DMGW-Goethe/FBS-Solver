import numpy as np
import matplotlib.path as mplPath # to check if a point is inside a polygon
import matplotlib.pyplot as plt	# for the contour plot
import pandas as pd

from . import plotting

# find the nearest point in the rho_c-phi_c diagram to a given input point
def findnearest_to_stabcurve(sol, stabcurve_point, index1, index2):

	nearestpointMR = [0,0]	# hold one point in the MR diagram
	distance = 1.0
	min_distance = 100.0

	# iterate through the whole array to find the smallest distance/ point nearest to stancurve_point
	for i in range(len(sol)):
		# compute the distance between a point in the sol array and the given stabcurve_point:
		distance = np.sqrt( (stabcurve_point[0]-sol[i][index1])**2 + (stabcurve_point[1]-sol[i][index2])**2 ) # set distance
		# find the nearest point:
		if (distance < min_distance):
			nearestpointMR = [sol[i][0], sol[i][4]]  # [M, R] point
			min_distance = distance  # update the current smallest distance

	# output the nearest point:
	return nearestpointMR


def refactorList(list, n):
	return [list[i:i+n] for i in range(0, len(list), n)]


def find_stabcurve_vectorfield(lines, df, indices, curve_index):

	#lines = sorted(lines, key=len, reverse=True)
	#stab_curve = lines[0]
	rhoVals = sorted(np.unique(df[:,indices["rho_0"]]))
	phiVals = sorted(np.unique(df[:,indices["phi_0"]]))
	numStarsRho = len(rhoVals)
	numStarsPhi = len(phiVals)
	E_0_max = max(phiVals)

	# filter out a pure boson star:
	pureBS = plotting.searchPureBS(df, indices)
	#pureBS = np.array([sol for sol in pureBS if np.abs(sol[indices[Y]]) < 10 * ylim[1] and np.abs(sol[indices[Y]]) > ylim[0]])
	pureBS_MT = pureBS[:,indices["M_T"]]
	pureBS_E0 = pureBS[:,indices["phi_0"]]
	
	M_T_index = 0
	for i in range(1,len(pureBS_MT)):
		M_T_index = i
		if pureBS_MT[i-1] > pureBS_MT[i]:
			break
		
	#print(M_T_index)
	#print(len(pureBS_MT))
	#print(pureBS_MT)
	if M_T_index >= (len(pureBS_MT)-2):
		# we have only a section of the stability curve (i.e. stab curve does NOT touch the y-axis)
		lines = sorted(lines, key=len, reverse=True)
		stab_curve = lines[curve_index]
		#stab_curve = np.append(np.array([[5.0, 5.1], [6.0, 6.1]]), stab_curve, axis=0)
		stab_curve = np.append(np.array([[0.0, E_0_max], [stab_curve[0][0], E_0_max]]), stab_curve, axis=0) # add line connecting the max E_0-value with rho_c>0 to the point with E
		return stab_curve
	else:	# the stability curve is either included completely or is cut into two pieces
		lines = sorted(lines, key=len, reverse=True)
		stab_curve = lines[curve_index]
		stab_curve_phivals = stab_curve[:,1]
		#print("stab_curve entries")
		#print(stab_curve)
		#print(stab_curve_phivals)
		if max(stab_curve_phivals) < 0.999*E_0_max:
			# the stab curve is completely included inside the stability plot
			lines = sorted(lines, key=len, reverse=True)
			stab_curve = lines[curve_index]
			return stab_curve
		else:
			# the stab curve is cut into two parts. Here we must stitch parts of the stabilitycurve together to one:
			# first find the stabcurve which touches the x-axis at the correct point:
			# filter out a pure boson star:
			pureNS = plotting.searchPureNS(df, indices)
			#pureBS = np.array([sol for sol in pureBS if np.abs(sol[indices[Y]]) < 10 * ylim[1] and np.abs(sol[indices[Y]]) > ylim[0]])
			pureNS_MT = pureNS[:,indices["M_T"]]
			pureNS_rho0 = pureNS[:,indices["rho_0"]]
			M_T_index2 = 0
			for i in range(1,len(pureNS_MT)):
				M_T_index2 = i
				if pureNS_MT[i-1] > pureNS_MT[i]:
					break
			
			last_stable_rhoc = pureNS_rho0[M_T_index2]
			# second, find the stabcurve which touches the y-axis at the right point
			last_stable_phic = pureBS_E0[M_T_index]
			# here we must stitch parts of the stabilitycurve together
			# loop through each stab line to find the correct one and then take the longest of the found ones:
			phic_axis_lines = []
			rhoc_axis_lines = []
			print(lines)
			for line in lines:
				if abs(line[0,1]/last_stable_phic -1.) < 0.15: #this line touches the y-axis at the right point
					phic_axis_lines.append(line)
				if abs(line[-1,0]/last_stable_rhoc -1.) < 0.15: #this line touches the y-axis at the right point
					rhoc_axis_lines.append(line)
			phic_axis_lines = sorted(phic_axis_lines, key=len, reverse=True)
			rhoc_axis_lines = sorted(rhoc_axis_lines, key=len, reverse=True)
			phi_c_line = phic_axis_lines[0]
			rho_c_line = rhoc_axis_lines[0]
			# stitch both lines together to get one stability line:
			return np.append(phi_c_line, rho_c_line, axis=0)




	# check if the whole stabcurve is included in the range:
	#points in stabcurve have with [rho_c, Phi_c]
	#for i in range(len(lines)):
	#	#lines[i][index_in_line][(rho or phi)]
	#	tmpline = lines[i]
	#	tmpline[:][1]

	lines = sorted(lines, key=len, reverse=True)
	stab_curve = lines[curve_index]
	return stab_curve


# Calculates the stability curve according to the paper by Henriques et al. "Stabiliy of Boson-Fermion Stars"
def calc_stability_curve(df, indices, debug = False, curve_index=0, rhoCutoff = 1, phiCutoff = 1, isFPS = False):
	# first, get some useful quantities from the list of all values
	rhoVals = sorted(np.unique(df[:,indices["rho_0"]]))
	phiVals = sorted(np.unique(df[:,indices["phi_0"]]))

	numStarsRho = len(rhoVals)
	numStarsPhi = len(phiVals)

	data_frame = pd.DataFrame(df[:, [indices['M_T'], indices['N_F'], indices['rho_0'], indices['phi_0']]], columns=['M_T', 'N_F', 'rho_0', 'phi_0'])
	#print(data_frame.pivot_table(values='M_T', index='rho_0', columns='phi_0'))
	factoredMasses = data_frame.pivot_table(values='M_T', columns='rho_0', index='phi_0')
	factoredFermionNumbers = data_frame.pivot_table(values='N_F', columns='rho_0', index='phi_0')

	# use the np.gradient function to calculate all gradient values at once
	massesGradient = np.gradient(factoredMasses, phiVals, rhoVals, edge_order=2)
	fermionNumbersGradient = np.gradient(factoredFermionNumbers, phiVals, rhoVals, edge_order=2)

	# prepare some lists that are then used to save relevant quantities
	deriv_array_2D = np.zeros(numStarsPhi * numStarsRho) # this is the derivative of the fermion number in the direction of constant mass
	massDerivs = np.zeros(numStarsPhi * numStarsRho) # abs of the mass gradient
	nfDerivs = np.zeros(numStarsPhi * numStarsRho) # abs of the fermion number gradient

	deriv_array_2D = refactorList(deriv_array_2D, numStarsRho)
	massDerivs = refactorList(massDerivs, numStarsRho)
	nfDerivs = refactorList(nfDerivs, numStarsRho)

	# loop over all values in order to fill in the above allocated lists
	for iRho in range(numStarsPhi):
		for iPhi in range(numStarsRho):
			grad_M = [0.0,0.0]
			grad_Nf = [0.0,0.0]

			grad_M[0], grad_M[1] = massesGradient[0][iRho][iPhi], massesGradient[1][iRho][iPhi]
			grad_Nf[0], grad_Nf[1] = fermionNumbersGradient[0][iRho][iPhi], fermionNumbersGradient[1][iRho][iPhi]

			# define a vector perpendicular to the derivative (gradiant) of M and normalize it:
			perp = np.array([grad_M[1], -grad_M[0]])
			perp /= np.linalg.norm(perp)

			# calc dot product of Nf gradient wit line perp to mass:
			# this is the directional derivative of Nf in the direction where M=const
			deriv_array_2D[iRho][iPhi] = perp.dot(grad_Nf)
			massDerivs[iRho][iPhi] = np.linalg.norm(grad_M)
			nfDerivs[iRho][iPhi] = np.linalg.norm(grad_Nf)


	# now that we have all values ready, we need to find the curve in deriv_array_2D that has zero value
	# also, below I am throwing out the values at the boundary (that's why I do phiVals[1:]) because the derivative is just not well behaved at the boundarys
	Yy, Xx = np.meshgrid(rhoVals[rhoCutoff:], phiVals[phiCutoff:])
	deriv_array_2D = [temp[rhoCutoff:] for temp in deriv_array_2D[phiCutoff:]]

	# plt.contour is perfect to find the curve with zero value
	fig = plt.figure(999)
	ax = fig.add_subplot(111)
	contours = ax.contour(Yy, Xx, deriv_array_2D, colors='black', levels=[0.00], antialiased = False)

	# extract all lines that were found with zero value
	lines = []
	for line in contours.collections[0].get_paths():
		lines.append(line.vertices)

	if isFPS:
		stab_curve = find_stabcurve_vectorfield(lines, df, indices, curve_index)
		#lines = sorted(lines, key=len, reverse=True)
		#stab_curve = lines[curve_index]
	else:
		lines = sorted(lines, key=len, reverse=True)
		stab_curve = lines[curve_index]

	#print(lines)
	#print(lines[curve_index])
	#stab_curve = np.append(np.array([[5.0, 5.1], [6.0, 6.1]]), stab_curve, axis=0)
	#print(stab_curve)
	#exit()
	ax.clear()
	plt.close(fig)

	# if debug = True, then also plot some of the intermediate quantities
	if(debug):
		Y, X = np.meshgrid(rhoVals, phiVals)

		plt.pcolormesh(Yy, Xx, (deriv_array_2D), cmap = "PRGn", alpha=0.8)
		plt.colorbar()
		plt.clim([-1,1])
		plt.contour(Yy, Xx, deriv_array_2D, colors='black', levels=[0.00], antialiased = False, linewidths = 4)
		plt.plot(stab_curve[:,0], stab_curve[:,1], linestyle="--", color = "red", linewidth = 4, label = "extracted stability line")
		plt.title("final deriv values")
		plt.legend(loc = "upper right")
		plt.show()

		plt.contour(Yy, Xx, deriv_array_2D, colors='black', levels=[0.00], antialiased = False, linewidths = 4)
		plt.plot(stab_curve[:,0], stab_curve[:,1], linestyle="--", linewidth = 4, color = "red", label = "extracted stability line")
		plt.pcolormesh(Y, X, (massDerivs), cmap = "viridis", alpha=0.8)
		plt.clim([0, 1000])
		plt.colorbar()
		plt.title("mass derivs")
		plt.legend(loc = "upper right")
		plt.show()

		plt.contour(Yy, Xx, deriv_array_2D, colors='black', levels=[0.00], antialiased = False, linewidths = 4)
		plt.plot(stab_curve[:,0], stab_curve[:,1], linestyle="--", linewidth = 4, color = "red", label = "extracted stability line")
		plt.pcolormesh(Y, X, (nfDerivs), cmap = "viridis", alpha=0.8)
		plt.clim([0, 1000])
		plt.colorbar()
		plt.title("nf derivs")
		plt.legend(loc = "upper right")
		plt.show()

	#exit()
	return stab_curve


# use the stability curve to filter out which FBS solutions are inside the stability region and which are not:
# return the stable configurations only:
def filter_stab_curve_data(sol_array, indices, stab_curve):

	# extend the stability curve to make it a closed polygon
	# Then, use it to search for all star configurations inside of the stability curve polygon
	# extend the stab_curve to a polygon
	# append additional elements:
	#stab_curve_polygon = np.append(stab_curve, np.array([[-0.5,-0.5]]), axis=0)  # chose a point in the lower left corner which is guaranteed to produce a polygon with the remaining curve
	stab_curve_polygon = np.append(stab_curve, np.array([[stab_curve[-1][0], 0.0]]), axis=0)
	stab_curve_polygon = np.append(stab_curve_polygon, np.array([[0.0, 0.0]]), axis=0)
	stab_curve_polygon = np.append(stab_curve_polygon, np.array([[0.0, stab_curve[0][1]]]), axis=0)

	#print(stab_curve_polygon)
	bbPath = mplPath.Path(stab_curve_polygon) # convert stab curve so that mplPath can use it

	# construct check if the stars are inside of the stability curve
	# (stab curve is a polygon and we check which points are inside it!)
	filtered_data = [] # create dummy output array

	# iterate through the array to see which point is inside the stab curve:
	for i in range(0,len(sol_array)):
		point = [sol_array[i][indices['rho_0']], sol_array[i][indices['phi_0']]]  # point with [rho_c, Phi_c]
		accuracy = 1e-10
		# check if point is inside:
		if (bbPath.contains_point(point,radius=accuracy) or bbPath.contains_point(point,radius=-accuracy)):
			tmplist = sol_array[i]
			filtered_data.append(tmplist) # append points inside the polygon (i.e. stable NS configurations) to the output array

	# return the filtered data which holds ONLY the stable configurations:
	# but before that, convert the python list to a numpy ndarray
	# (we need to do this to ensure compatibility with other plotting functions)
	return np.array(filtered_data)


# converts a stability curve in rho_c-phi_c space into MR space, so that a polygon can be constructed from it
def stab_curve_to_MR(sol_array, indices, old_stabcurve, num_rho_stars, num_phi_stars):

	MR_stabcurve = [] # list which holds arrays with 2 entries

	# iterate through the stability curve in the rhoc-phi_c diagram and find the nearest point to each stab curve segment respectively:
	for i in range(len(old_stabcurve)):
		MRpoint = findnearest_to_stabcurve(sol_array, old_stabcurve[i], indices['rho_0'], indices['phi_0'])  # [M, R] point
		MR_stabcurve.append(MRpoint) # append the nearest point

	# now the lines of rho_c=0 and phi_c=0 need to be added to the curve:
	# line where rho_c = 0:
	MR_stabcurve_rho0 = []
	for j in range(num_phi_stars):
		index = j*num_rho_stars
		MRpoint = [sol_array[index][indices['M_T']], sol_array[index][indices['R_F_0']]] # [M, R] point
		# iterate until the point is not on the stability curve anymore
		if (sol_array[index][indices['phi_0']] > old_stabcurve[len(old_stabcurve)-1][1]):   #check for phi_c < phi_c stabcurve
			break

		MR_stabcurve_rho0.append(MRpoint)

	# line where phi_c = 0:
	MR_stabcurve_phi0 = []
	for j in range(num_rho_stars):
		index = j
		MRpoint = [sol_array[index][indices['M_T']], sol_array[index][indices['R_F_0']]] # [M, R] point
		if (sol_array[index][indices['rho_0']] > old_stabcurve[0][0]):   #check for rho_c < rho_c stabcurve
			break

		MR_stabcurve_phi0.append(MRpoint)

	# append the MR curves as phi_c=0 and rho_c=0 to the stability curve:
	tmp1 = MR_stabcurve_rho0[::-1] # reverse element order to get the correct item order for appending later
	for k in range(len(MR_stabcurve_rho0)):
		MR_stabcurve.append(tmp1[k])

	for j in range(len(MR_stabcurve_phi0)):
		MR_stabcurve.append(MR_stabcurve_phi0[j])

	return np.array(MR_stabcurve) # return the stability curve in the MR diagram
