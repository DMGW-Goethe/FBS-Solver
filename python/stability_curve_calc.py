import numpy as np
import matplotlib.pyplot as plt	# for the contour plot


def create_3D_array(num_i, num_j, len_inner_arr):
	tmp = []
	for i in range(num_i):
		tmp2 = []
		for j in range(num_j):
			tmp2.append([0.]*len_inner_arr)
		tmp.append(tmp2)
	return tmp


# create a value array with padding ghost cells
def create_ghost_cell_array(sol_array, num_rho_stars, num_phi_stars, stencil_order):

	# create the array which holds data AND ghost cells:
	data_array_w_ghost_cells = create_3D_array(num_rho_stars+2*stencil_order, num_phi_stars+2*stencil_order, len(sol_array[0]))

	# extract the Nb and Nf values and arrange them into a 2D array:
	# in each entry we now hold all values corresponding to one star
	for i in range(num_rho_stars):
		for j in range(num_phi_stars):
			for k in range(len(sol_array[0]) ):
				tmp = sol_array[i+ num_rho_stars*j][k]
				data_array_w_ghost_cells[i+stencil_order][j+stencil_order][k] = tmp
				#print(data_array_w_ghost_cells[i][j])
	
	# fill the corners of the ghost cell-domain:
	for i in range(stencil_order):
		for j in range(stencil_order):
			# lower left corner (use the upper left corner value):
			data_array_w_ghost_cells[i][j] = np.array(sol_array[0])
			# upper right corner:
			data_array_w_ghost_cells[i+stencil_order+num_rho_stars][j+stencil_order+num_phi_stars] = np.array(sol_array[num_rho_stars*num_phi_stars-1])
			# lower right corner:
			data_array_w_ghost_cells[i+stencil_order+num_rho_stars][j] = np.array(sol_array[num_rho_stars-1])
			# upper left corner:
			data_array_w_ghost_cells[i][j+stencil_order+num_phi_stars] = np.array(sol_array[num_rho_stars*num_phi_stars-num_rho_stars])

	# set the sides without the corners:
	# left and right side:
	for i in range(stencil_order):
		for j in range(num_phi_stars):
			# left side:
			data_array_w_ghost_cells[i][j+stencil_order] = sol_array[num_rho_stars*j]
			# right side:
			data_array_w_ghost_cells[i+stencil_order+num_rho_stars][j+stencil_order] = sol_array[num_rho_stars*j+num_rho_stars-1]
	
	# upper and lower side:
	for i in range(num_rho_stars):
		for j in range(stencil_order):
			# lower side:
			data_array_w_ghost_cells[i+stencil_order][j] = sol_array[i]
			# upper side:
			data_array_w_ghost_cells[i+stencil_order][j+stencil_order+num_phi_stars] = sol_array[num_rho_stars*num_phi_stars-num_rho_stars+i]

	return data_array_w_ghost_cells


# central (symmetric) stencil
def central_stencil(vals, i, j, order, val_index):

	stencil2nd = [-0.5, 0.0, 0.5]
	stencil4th = [1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0]
	stencil6th = [-1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0]

	stencil = []
	if order == 2:
		stencil = stencil2nd
	elif order == 4:
		stencil = stencil4th
	elif order == 6:
		stencil = stencil6th
	else:
		print("CAUTION: wrong stencil order")

	derivIdir = 0
	derivJdir = 0
	for iTemp in range(-order//2, order//2 + 1):
		derivIdir += vals[i + iTemp][j][val_index] * stencil[iTemp + order//2] # derivative in rho_c dir
		derivJdir += vals[i][j + iTemp][val_index] * stencil[iTemp + order//2] # derivative in phi_c dir
    
	delta_rho = vals[i + 1][j][1] - vals[i][j][1]
	delta_phi = vals[i][j + 1][2] - vals[i][j][2]

	#print(delta_rho, delta_phi)
	derivIdir /= delta_rho	# derivative in rho_c dir
	derivJdir /= delta_phi	# derivative in phi_c dir
	return derivIdir, derivJdir # rho-dir, phi-dir

# forward stenctils
def forward_stencil(vals, i, j, order, val_index):
    stencil1st = [-1.0, 1.0]
    stencil2nd = [-3.0/2.0, 2.0, -1.0/2.0]
    stencil3rd = [-11.0/6.0, 3.0, -3.0/2.0, 1.0/3.0]
    stencil4th = [-25.0/12.0, 4.0, -3.0, 4.0/3.0, -1.0/4.0]

    stencil = []
    if order == 1:
        stencil = stencil1st
    elif order == 2:
        stencil = stencil2nd
    elif order == 3:
        stencil = stencil3rd
    elif order == 4:
        stencil = stencil4th

    derivIdir = 0
    derivJdir = 0
    for iTemp in range(0, order + 1):
        # print(iTemp)
        derivIdir += vals[i + iTemp][j][val_index] * stencil[iTemp]
        derivJdir += vals[i][j + iTemp][val_index] * stencil[iTemp]

    derivIdir /= vals[i + 1][j][1] - vals[i][j][1]
    derivJdir /= vals[i][j + 1][2] - vals[i][j][2]
    return derivIdir, derivJdir


# calculate the stability curve using the definition in arXiv:2006.08583v2
# by computing partial derivatives of Nb and Nf
# data must be parametrized in terms of rho_c and phi_c!
def calc_stability_curve(sol_array, num_rho_stars, num_phi_stars, stencil_order):

	# create a 2D array to hold the values of the partial derivatives which are used to find the stability curve:
	deriv_array_2D = create_3D_array(num_rho_stars, num_phi_stars, 1) # 3D array with the last index having length=1 is just a 2D array

	# extract the Nb and Nf values (actually all physical values) and arrange them into a 2D array:
	data_array_2D_w_ghost_cells = create_ghost_cell_array(sol_array, num_rho_stars, num_phi_stars, stencil_order)

	# now use this array filled with the wanted values AND padding ghost cells togehter with a stencil of wanted order:
	for i in range(num_rho_stars):
		for j in range(num_phi_stars):
			irho = i + stencil_order
			jphi = j + stencil_order
			# holds the partial derivative (gradient):
			grad_M = [0.0,0.0]
			grad_Nf = [0.0,0.0]
			
			# (last function argument: 0=totalmass, 5=fermion number, 7=boson number)
			grad_M[0], grad_M[1] = central_stencil(data_array_2D_w_ghost_cells, irho, jphi, stencil_order, 0) #calc derivative of M in rho and phi dir
			grad_Nf[0], grad_Nf[1] = central_stencil(data_array_2D_w_ghost_cells, irho, jphi, stencil_order, 5) #calc deriv. of Nf in rho and phi dir

			# define a vector perpendicular to the derivative (gradiant) of M and normalize it:
			perp = np.array([grad_M[1], -grad_M[0]])
			perp /= np.linalg.norm(perp)

			# calc dot product of Nf gradient wit line perp to mass:
			# this is the directional derivative of Nf in the direction where M=const
			deriv_array_2D[j][i] = perp.dot(grad_Nf)

	# obtain all points where deriv_array is (almost) zero/ has a minimum and add these points to the stability curve:
	stab_curve = [] # init empty array to hold the stability line

	# obtain the contour line using the plt.contour function!:
	rhogrid = []
	phigrid = []
	# create a grid in rho_c and phi_c using the available data:
	for i in range(num_rho_stars):
		rhogrid.append(data_array_2D_w_ghost_cells[i+stencil_order][0][1])
	for j in range(num_phi_stars):
		phigrid.append(data_array_2D_w_ghost_cells[0][j+stencil_order][2])
	# now call the function:
	Y, X = np.meshgrid(rhogrid, phigrid)
	contours = plt.contour(Y, X, deriv_array_2D, colors='black', levels=[0.00], antialiased = True)

	#plt.imshow(deriv_array_2D, extent=[0.0, 0.008, 0.0, 0.14], origin='lower', cmap='turbo', aspect='auto', interpolation='none', alpha=0.8)

	# extract the contour lines:
	lines = []
	for line in contours.collections[0].get_paths():
		lines.append(line.vertices)

	# select only the longest line (this is the one that we want):
	longestindex = 0
	maxlen = 0
	for k in range(len(lines)):
		print(len(lines[k]))
		if len(lines[k]) > maxlen:
			maxlen = len(lines[k])
			longestindex = k

	# the wanted stability curve is the line which has the longest length (in case there would be multiple lines)
	stab_curve = lines[longestindex]

	plt.clf()	# to prevent the contour to be plotted (also we need this to remove the small "islands" in the plot)
	# return the (raw/unprocessed) stability curve:
	return stab_curve


# use the stability curve to filte out which FBS solutions are inside the stability region and which are not:
# return the stable configurations only:
def filter_stab_curve_data():


	# call stability curve finder:


	# construct check if the stars are inside of the stability curve
	# (stab curve is a polygon and we check which points are inside it!)

	filtered_data = []

	# return the filtered data which holds ONLY the stable configurations
	return filtered_data
