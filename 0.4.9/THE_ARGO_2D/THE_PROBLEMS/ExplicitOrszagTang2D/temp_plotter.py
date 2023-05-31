## temp_plotter_2d.py
## by Ryan Farber 21 January  2015
## Last modified: 01 May      2015
"""
The purpose of this program is to load fluid variables processed by
The Argo 2D from a saved state file and to plot the 2D variable.

Indexing of "the_data" depends on what dimension of "save_state"
was used; the below disambiguates for data from The Argo 2D.

the_data[ 0] = cycles
the_data[ 1] = u       # x-component of velocity
the_data[ 2] = v       # y-component of velocity
the_data[ 3] = cn      # concentration (tracer density)
the_data[ 4] = p       # pressure
the_data[ 5] = src     # src term of pressure poisson eqn
the_data[ 6] = Fx      # x-comp of external force
the_data[ 7] = Fy      # y-comp of external force
the_data[ 8] = Bx      # x-comp of magnetic field
the_data[ 9] = By      # y-comp of magnetic field
the_data[10] = A       # vector potential
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
plot_file_ext = ".tiff"

os.chdir(cur_path + "/StateFiles")
datafiles = glob.glob("*cycle*" + data_file_ext)

for datafile in datafiles:
    the_data = cPickle.load(open(datafile, "rb"))
    os.chdir(cur_path + "/Plots")
    
    plt_name1  = "vMag" + datafile[:-len(data_file_ext)] + plot_file_ext
    plt_name2  = "bEng" + datafile[:-len(data_file_ext)] + plot_file_ext
    
    X = np.linspace(my_fluid.XMIN, my_fluid.XMAX, my_fluid.NX)
    Y = np.linspace(my_fluid.YMIN, my_fluid.YMAX, my_fluid.NY)
    Y,X = np.meshgrid(X,Y)
    X = X[1:-1,1:-1]
    Y = Y[1:-1,1:-1]

    u  = the_data[1][1:-1,1:-1]
    v  = the_data[2][1:-1,1:-1]
    Bx = the_data[8][1:-1,1:-1]
    By = the_data[9][1:-1,1:-1]
                
    vMag = np.sqrt(u**2 + v**2)
    bEng = 0.5*(Bx**2 + By**2)

    plt_title1 = "Magnitude of Velocity Colormesh at Cycle " + str(the_data[0])
    plt_title2 = "Magnetic Field Energy Density Colormesh at Cycle " + str(the_data[0])
            
    plt.pcolormesh(X,Y, vMag)
    plt.colorbar()
                            
    plt.xlabel("x"); plt.ylabel("y")
    plt.title(plt_title1)
    plt.savefig(plt_name1)

    plt.clf()

    plt.pcolormesh(X,Y, bEng)
    plt.colorbar()
                            
    plt.xlabel("x"); plt.ylabel("y")
    plt.title(plt_title2)
    plt.savefig(plt_name2)

    plt.clf()
    
    os.chdir(cur_path + "/StateFiles")
# end for
        
## end temp_plotter_2d.py
