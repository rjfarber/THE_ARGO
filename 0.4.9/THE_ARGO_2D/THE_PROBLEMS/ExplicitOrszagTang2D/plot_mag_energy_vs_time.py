## plot_mag_energy_vs_time.py
## by Ryan Farber 28 April 2015
## Last modified: 28 April 2015
"""
The purpose of this program is to load fluid variables processed by
The Argo 2D from a saved state file and to plot the 2D variable.

Indexing of "the_data" depends on what dimension of "save_state"
was used; the below disambiguates for data from The Argo 2D.

the_data[0] = cycles
the_data[1] = u       # x-component of velocity
the_data[2] = v       # y-component of velocity
the_data[3] = cn      # concentration (tracer density)
the_data[4] = p       # pressure
the_data[5] = src     # src term of pressure poisson eqn
the_data[6] = Fx      # x-comp of external force
the_data[7] = Fy      # y-comp of external force
the_data[8] = Bx      # x-comp of magnetic field
the_data[9] = By      # y-comp of magnetic field
"""
import glob
import numpy as np
import matplotlib.pyplot as plt; fig = plt.figure()
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import sys; sys.path.insert(0,"../../..")
from the_fluid import The_Fluid
import os; cur_path = os.getcwd(); os.chdir(cur_path + "/StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Parameters
data_file_ext = ".p"
plot_file_ext = ".png"

os.chdir(cur_path + "/StateFiles")
datafiles = glob.glob("*cycle*" + data_file_ext)

cycles = []
mag_energy = []
for datafile in datafiles:
    the_data = cPickle.load(open(datafile, "rb"))
    
    cycles.append(the_data[0])
    
    Bx = the_data[8]
    By = the_data[9]

    Bx = Bx[1:-1,1:-1]
    By = By[1:-1,1:-1]

    mag_energy.append(sum(sum(Bx**2 + By**2)))

    if max(cycles) > 200:
        break
# end for
os.chdir(cur_path + "/Plots")
plt.plot(cycles, mag_energy)
plt.savefig("magnetic_energy_over_time.png")
        
## end plot_mag_energy_vs_time.py
