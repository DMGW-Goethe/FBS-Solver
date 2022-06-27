import numpy as np

# loads the custum FBS data values (but only the numbers and not the letters):
def load_MRPhi_data(myfilename):
	data = np.loadtxt(myfilename) #, delimeter = ' ')
	return data
