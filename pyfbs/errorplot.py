from . import stability_curve as scc # import custom python file which computes the stability curve
import numpy as np
# a file which contains the functions to comute the relative errors of two quantities. To be used by plot_error_comparison_effsys_fullsys()

# function to compute min, max and average relative error for one configuration:
def calc_relative_err_min_max_median_sigma(data1, data2):
    
	# filter out values in the data1, data2 arrays which are obviously wrong e.g. NaN or very small/large values:
#	data1_filtered = np.zeros(0)
#	data2_filtered = np.zeros(0)
#	for i in range(len(data1)):
#		if ( np.isnan(data1[i]) or np.isnan(data2[i]) or (data1[i] < 1e-10) or (data2[i] < 1e-10) or (data1[i] > 1e10) or (data2[i] > 1e10)):
#			print(data1[i], data2[i])
#		else:
#			np.append(data1_filtered,data1[i])# append values which are not broken
#			np.append(data2_filtered,data2[i])
#	print(data1_filtered)
#	error_list = abs(data1_filtered-data2_filtered)	# absolute error abs(data2)
#	error_list = error_list / abs(data1_filtered)	# relative error
	
	error_list = abs(abs(data1)-abs(data2))	# absolute error abs(data2)
	error_list = error_list / abs(data1)	# relative error
	#max_err = max(error_list)
	#max_err = min(max_err, 1.)	# failsave to
	min_err = min(error_list)
	# first we have to sort the list to find the median value and the 1 and 2 sigma deviations:
	# median 50% of values are smaller, 1sigma: 68% is smaller, 2sigma: 95% is smaller
	error_list.sort()
	max_err = error_list[-2] # for now exclude the largest value by hand since it likely is NaN or ridiculously large
	median_err = error_list[(int)(50./100 * len(error_list))]
	one_sigma  = error_list[(int)(68./100 * len(error_list))]
	two_sigma  = error_list[(int)(95./100 * len(error_list))]
	#print(error_list)

	return min_err, max_err, median_err, one_sigma, two_sigma


# currently not in use!
# filter out solutions with radii > 100 km (these objects are not Neutron Stars anyways) and using other critria (e.g. NaNs)
def filter_errorplot_data(data1, data2, indices1, indices2):
	newd1 = []
	newd2 = []
	# now only apend values to the arrays if they do not match one of the filter criteria
	for i in range(len(data1)): # data1 and data2 should have the same length
		filter_out_condition = False
		filter_out_condition = (data1[i][indices1["R_F"]] > 50.) or (data2[i][indices2["R_F_0"]] > 50.)

		if not filter_out_condition:
			newd1.append(data1[i])
			newd2.append(data2[i])

	return np.array(newd1), np.array(newd2)


def calc_error_curves(data_in, indices_fullsys, indices_effsys, myindex, filter_stars=False):

	if ((len(data_in) % 2) == 1):	# include a check because there must be pairs of data: full system + effective system
		print("need an even amount of data files!\n With even files containing data of the full system and uneven from the effective system!")
	
	# pre-processing of data files:
	# filter the data
	myworkdata = [None]*len(data_in)
	stab_curve = [None]*len(data_in)
	# compute all stability curves using the solution from the full system:
	for k in range(0,len(data_in),2):
		stab_curve[k] = scc.calc_stability_curve(data_in[k], indices_fullsys, debug = False, curve_index = 0)
		myworkdata[k] = scc.filter_stab_curve_data(data_in[k], indices_fullsys, stab_curve[k]) # filter the stars inside the stability region
		myworkdata[k+1] = scc.filter_stab_curve_data(data_in[k+1], indices_fullsys, stab_curve[k]) # filter the stars inside the stability region
	
	# filter out data using other criteria:
	if (filter_stars):
		for l in range(0,len(data_in),2):
			myworkdata[l], myworkdata[l+1] = filter_errorplot_data(myworkdata[l], myworkdata[l+1], indices_fullsys, indices_effsys)
        
	# compute the reative min/max/average errors
	minerr = []	# def helper arrays:
	maxerr = []
	medianerr = []
	one_sigmaerr = []
	two_sigmaerr = []

	# calculate the errors of a given quantity (given by the index):
	for j in range(0,len(myworkdata),2):
		data1 = myworkdata[j][1:,indices_fullsys[myindex]]	# tmp variables
		data2 = myworkdata[j+1][1:,indices_effsys[myindex]]
		# here, the errors are calculated:
		minerr_tmp, maxerr_tmp, medianerr_tmp, one_sigmaerr_tmp, two_sigmaerr_tmp = calc_relative_err_min_max_median_sigma(data1, data2)
		minerr.append(minerr_tmp)
		maxerr.append(maxerr_tmp)
		medianerr.append(medianerr_tmp)
		one_sigmaerr.append(one_sigmaerr_tmp)
		two_sigmaerr.append(two_sigmaerr_tmp)


	return minerr, maxerr, medianerr, one_sigmaerr, two_sigmaerr