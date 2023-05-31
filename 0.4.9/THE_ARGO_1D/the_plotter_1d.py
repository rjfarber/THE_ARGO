## the_plotter_1d.py
## by Ryan Farber 21 January 2015
## Last modified: 26 January 2015
"""
The purpose of this program is to load fluid variables processed by
The Argo 1D from a saved state file and to plot the 1D variable.

  NOTE: This file must be copied by a bash script to the specific
Problem file folder it shall be used in to work.

Indexing of "the_data" depends on what dimension of "save_state"
was used; the below disambiguates for data from The Argo 1D.

the_data[0] = cycles
the_data[1] = u or T
"""
import glob
import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path.insert(0,"../../..")
from the_fluid import The_Fluid
import os; cur_path = os.getcwd(); os.chdir(cur_path + "/StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Parameters
data_ext = my_fluid.DATA_FILE_EXT
plot_ext = my_fluid.PLOT_FILE_EXT

pfiles = glob.glob("*cycle*.p")
for pfile in pfiles:
    the_data = cPickle.load(open(pfile, "rb"))

    plt_title = "Cycle " + str(the_data[0])
    plt_name  = pfile[:-len(data_file_ext)] + plot_file_ext
    if "my_fluid.PLOT_TYP" in vars() or "my_fluid.PLOT_TYP" in globals():
        ylbl = my_fluid.PLT_TYP
    else:
        ylbl = "velocity"
    
    X = np.linspace(my_fluid.XMIN, my_fluid.XMAX, my_fluid.NX)
    
    os.chdir(cur_path + "/Plots")
    plt.plot(X, the_data[1])
    plt.xlim(my_fluid.XMIN, my_fluid.XMAX); plt.ylim(-2.5, 2.5)
    plt.xlabel('x'); plt.ylabel(ylbl)
    plt.title(plt_title)
    plt.savefig(plt_name)

    plt.clf()
    os.chdir(cur_path + "/StateFiles")
# end for
        
## end the_plotter_1d.py
