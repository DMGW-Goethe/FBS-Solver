import numpy as np
import matplotlib.pyplot as plt	# for the contour plot

# create a value array with padding ghost cells
def create_ghost_cell_array(sol_array, num_rho_stars, num_phi_stars, stencil_order):

	#data_array_w_ghost_cells = np.zeros(shape=(num_rho_stars+2*stencil_order, num_phi_stars+2*stencil_order, len(sol_array[0])))
	#data_array_w_ghost_cells = [[[None]*len(sol_array[0])]*(num_phi_stars+2*stencil_order)]*(num_rho_stars+2*stencil_order)
	#data_array_w_ghost_cells = [[[0]*len(sol_array[0]) ]* (num_phi_stars+2*stencil_order)]* (num_rho_stars+2*stencil_order)

	# create the array which holds data AND ghost cells:
	test = []
	for i in range(num_rho_stars+2*stencil_order):
		test2 = []
		for j in range(num_phi_stars+2*stencil_order):
			test2.append([0.]*len(sol_array[0]))
		test.append(test2)

	#print(data_array_w_ghost_cells)
	data_array_w_ghost_cells = test
	#print(data_array_w_ghost_cells)

	#print(test)
	#test[0][0][1] = 1.
	#print(test)

	# extract the Nb and Nf values and arrange them into a 2D array:
	# in each entry we now hold all values corresponding to one star
	for i in range(num_rho_stars):
		for j in range(num_phi_stars):
			for k in range(len(sol_array[0]) ):
				tmp = sol_array[i+ num_rho_stars*j][k]
				data_array_w_ghost_cells[i+stencil_order][j+stencil_order][k] = tmp
				#print(data_array_w_ghost_cells[i][j])

	# prepare the array including ghost cells:
	# set the corners:
	# for i in range(num_rho_stars+2*stencil_order):
	# 	for j in range(num_phi_stars+2*stencil_order):
	# 		#print(i,j)
	# 		#print(sol_array[i+j*num_rho_stars][0], end=" ")
	# 		print(data_array_w_ghost_cells[i][j][0], end = " ")
	# 	print()
	# exit()

	
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

	#exit()
	#print(data_array_w_ghost_cells[:,:,0])

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


	# for i in range(num_rho_stars+2*stencil_order):
	# 	for j in range(num_phi_stars+2*stencil_order):
	# 		#print(i,j)
	# 		print(data_array_w_ghost_cells[i][j][0], end=" ")
	# 	print()

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

	# arrays to hold the Nb and Nf values of all points on the stability curve:
	Rpos = []
	Mpos = []	
	#deriv_array_2D = [[0]*(num_phi_stars)]*(num_rho_stars)

	test = []
	for i in range(num_rho_stars):
		test2 = []
		for j in range(num_phi_stars):
			test2.append([0.])
		test.append(test2)
	
	deriv_array_2D = test

	# extract the Nb and Nf values and arrange them into a 2D array:
	data_array_2D_w_ghost_cells = create_ghost_cell_array(sol_array, num_rho_stars, num_phi_stars, stencil_order)

	#for i in range(num_rho_stars):
	#	for j in range(num_phi_stars):
	#		
	#		print(data_array_2D_w_ghost_cells[i][j][1])
	#print(sol_array)
	#exit()

	# now use this array filled with the wanted values AND padding ghost cells togehter with a stencil of wanted order:
	
	for i in range(num_rho_stars):
		for j in range(num_phi_stars):
			irho = i + stencil_order
			jphi = j + stencil_order
			# holds the derivative
			grad_M = [0,0]
			grad_Nf = [0,0]
			
			# (last function argument: 0=totalmass, 5=fermion number, 7=boson number)
			grad_M[0], grad_M[1] = central_stencil(data_array_2D_w_ghost_cells, irho, jphi, stencil_order, 0) #calc derivative of M in rho and phi dir
			grad_Nf[0], grad_Nf[1] = central_stencil(data_array_2D_w_ghost_cells, irho, jphi, stencil_order, 7) #calc deriv. of Nf in rho and ohi dir
			
			#grad_M[0], grad_M[1] = forward_stencil(data_array_2D_w_ghost_cells, irho, jphi, stencil_order, 0) #calc derivative of M in rho and phi dir
			#grad_Nf[0], grad_Nf[1] = forward_stencil(data_array_2D_w_ghost_cells, irho, jphi, stencil_order, 5) #calc deriv. of Nf in rho and ohi dir

			# define a vector perpendicular to the derivative in M
			perp = np.array([grad_M[1], -grad_M[0]])
			perp /= np.linalg.norm(perp)

			# calc dot product of Nf gradient wit line perp to mass:
			deriv_array_2D[j][i] = perp.dot(grad_Nf)

	# obtain all points where deriv_array is (almost) zero/ has a minimum and add these points to the stability curve:
	# save the stability curve in terms of rho_c and Phi_c AND in terms of the corresponding M and R:
	stab_curve = []
	# deprecated! needs to be removed later
	for i in range(num_rho_stars-1):
		for j in range(num_phi_stars-1):

			if ((deriv_array_2D[i][j] * deriv_array_2D[i][j+1] < 0) or (deriv_array_2D[i][j] * deriv_array_2D[i+1][j] < 0)):
				# add this point to the stability curve:
				Rpos = data_array_2D_w_ghost_cells[stencil_order+i][stencil_order+j][4]	# Rvalue
				Mpos = data_array_2D_w_ghost_cells[stencil_order+i][stencil_order+j][0]	# Mvalue
				Rho_c_pos = data_array_2D_w_ghost_cells[stencil_order+i][stencil_order+j][1] # Rhoc value
				Phi_c_pos = data_array_2D_w_ghost_cells[stencil_order+i][stencil_order+j][2] # Phic value
				stab_curve.append([Rpos,Mpos, Rho_c_pos, Phi_c_pos])
	
	#print(deriv_array_2D)


	# TODO: move all the plotting related parts of the code to other functions, once the stability curve calculation is working correctly
	rhogrid = []
	phigrid = []

	for i in range(num_rho_stars):
		rhogrid.append(data_array_2D_w_ghost_cells[i+stencil_order][0][1])
	for j in range(num_phi_stars):
		phigrid.append(data_array_2D_w_ghost_cells[0][j+stencil_order][2])

	# obtain the contour line in a different way:

	
	Y, X = np.meshgrid(rhogrid, phigrid)
	contours = plt.contour(Y, X, deriv_array_2D, colors='black', levels=[0.00], antialiased = True)
	#plt.show()
	#plt.imshow(deriv_array_2D, extent=[0.0, 0.008, 0.0, 0.14], origin='lower', cmap='turbo', aspect='auto', interpolation='none', alpha=0.8)
	#plt.show()

	plt.clf()	# to prevent the contour to be plotted (i.e we need this to remove the small "islands" in the plot)

	M_array = test
	for i in range(num_rho_stars):
		for j in range(num_phi_stars):
			M_array[j][i] = data_array_2D_w_ghost_cells[i+stencil_order][j+stencil_order][0]

	contours2 = plt.contour(rhogrid, phigrid, M_array, colors=['purple', 'brown', 'red', 'green', 'orange'], levels=[0.6250, 1.0, 1.2, 1.4, 1.6])

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

	print(lines[longestindex][:,0])

	plt.plot(lines[longestindex][:,0], lines[longestindex][:,1], c='black', linewidth=2)

	plt.imshow(M_array, extent=[0.0, 0.008, 0.0, 0.14], origin='lower', cmap='turbo', aspect='auto', interpolation='none', alpha=0.8)
	plt.show()

	return stab_curve



def filter_stab_curve_data():


	# call stability curve finder:


	# construct check if the stars are inside of the stability curve
	# (stab curve is a polygon and we check which points are inside it!)

	filtered_data = []

	# return the filtered data which holds ONLY the stable configurations
	return filtered_data
