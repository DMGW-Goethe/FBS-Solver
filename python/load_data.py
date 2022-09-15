import numpy as np

# loads the custum FBS data values (but only the numbers and not the letters):
def load_MRPhi_data(myfilename):
    f = open(myfilename)
    l0 = next(f)
    labels = l0.replace('#', '').strip().split('\t')
    indices = dict([(labels[i].strip(), i) for i in range(len(labels))])
    f.close()
    data = np.loadtxt(myfilename) #, delimeter = ' ')
    return data, indices
