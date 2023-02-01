import load_data
import plotting_functions as pfuncts
import stability_curve_calc as scc
import numpy as np
import matplotlib.pyplot as plt

import importlib
importlib.reload(pfuncts)
importlib.reload(scc)

from matplotlib import rcParams
from matplotlib import rcParamsDefault
from matplotlib import rc

rcParams.update(rcParamsDefault)

rc('font', **{'family': 'CMU Sans Serif'})
#rc('font', **{'family': 'CMU Serif'})
#rcParams["mathtext.fontset"] = "cm"

# load in data:
importlib.reload(scc)
importlib.reload(pfuncts)

# load in all needed data files for different Lambda but const phi_c
filename1 = "../plots/TLN-line_phic-002_effectivesys-mu_0.500000_2513.274123.txt"
filename2 = "../plots/TLN-line_phic-002_fullsys-mu_0.500000_2513.274123.txt"

#for i in range(len(df)):
df1, indices = load_data.load_MRPhi_data(filename1)
df2, indices2 = load_data.load_MRPhi_data(filename2)

# here, define the plotting routines to make a nice looking plot:

# create/filter the data to be plotted:

filteredData1 = df1
allMasses1 = filteredData1[:,indices['M_T']]
allTlns1 = filteredData1[:,indices["lambda_tidal"]]
allDimlessTlns1 = allTlns1 / allMasses1**5
#--
filteredData2 = df2
allMasses2 = filteredData2[:,indices['M_T']]
allTlns2 = filteredData2[:,indices["lambda_tidal"]]
allDimlessTlns2 = allTlns2 / allMasses2**5

# ----------------------------------------------------------

# define plot configurations etc:
#ylim = [1, 1e6], xlim = [0, 3]
tickFontSize = 13

fig = plt.figure()
plt.xlim([1., 10.])
plt.ylim([1., 1e7])
plt.yscale("log")

plt.ylabel(r"Tidal Deformability $\Lambda$", fontsize = tickFontSize)
plt.xlabel(r"Total Gravitational Mass [M$_\odot$]", fontsize = tickFontSize)
#plt.tick_params(axis='both', which='major', labelsize=tickFontSize)
plt.grid(alpha=0.2, linestyle="-")

#for j in range(len(df)):
#plt.plot(allMasses[0], allDimlessTlns[0], label = "$\phi_c = 0.02$")
#plt.plot(allMasses[1], allDimlessTlns[1], label = "$\phi_c = 0.02$")

plt.plot(allMasses1, allDimlessTlns1, label = "$1 \phi_c = 0.02$")
plt.plot(allMasses2, allDimlessTlns2, label = "$2 \phi_c = 0.02$")

plt.show()
