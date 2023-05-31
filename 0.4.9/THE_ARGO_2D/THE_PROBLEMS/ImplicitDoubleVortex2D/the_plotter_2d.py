## the_plotter_2d.py
## by Ryan Farber 21 January 2015
## Last modified: 27 January 2015
"""
The purpose of this program is to load fluid variables processed by
The Argo 2D from a saved state file and to plot the 2D variable.

  NOTE: This file must be copied by a bash script to the specific
Problem file folder it shall be used in to work.

Indexing of "the_data" depends on what dimension of "save_state"
was used; the below disambiguates for data from The Argo 2D.

the_data[0] = cycles
the_data[1] = u       # x-component of velocity
the_data[2] = v       # y-component of velocity
the_data[3] = rho     # tracer density
the_data[4] = p       # pressure
the_data[5] = src     # src term of pressure poisson eqn
the_data[6] = F       # external force

Note that my_fluid.PLT_TYP must be constructed in the INPUTS* file
with care for proper operation of this file. Below are suggestions
(note: order of the keywords does not matter and spacing between
keywords similarly does not matter):

my_fluid.PLT_TYP = "quiver"               # quiver only plot of velocity
my_fluid.PLT_TYP = "quiver"+"pcolormesh"  # quiver plot of velocity on a
                                          # tracer density pcolormesh
my_fluid.PLT_TYP = "quiver"+"contour"     # quiver plot of velocity on
                                          # pressure contours
my_fluid.PLT_TYP = "surface"+"velocity"   # surface   plot of velocity
my_fluid.PLT_TYP = "surface"+"pressure"   # surface   plot of pressure
my_fluid.PLT_TYP = "wireframe"+"velocity" # wireframe plot of velocity
my_fluid.PLT_TYP = "wireframe"+"pressure" # wireframe plot of pressure
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

os.chdir(cur_path + "/Plots")
plotfiles = glob.glob("*cycle*" + plot_file_ext)

os.chdir(cur_path + "/StateFiles")
datafiles = glob.glob("*cycle*" + data_file_ext)

pfiles = glob.glob("*cycle*.p")
for datafile in datafiles:
    if (datafile[:-len(data_file_ext)] + plot_file_ext) in plotfiles:
        continue
    the_data = cPickle.load(open(datafile, "rb"))
    os.chdir(cur_path + "/Plots")

    plt_title = "Cycle " + str(the_data[0])
    plt_name  = datafile[:-len(data_file_ext)] + plot_file_ext
    
    X = np.linspace(my_fluid.XMIN, my_fluid.XMAX, my_fluid.NX)
    Y = np.linspace(my_fluid.YMIN, my_fluid.YMAX, my_fluid.NY)
    Y,X = np.meshgrid(X,Y)
    
    if   "quiver" in my_fluid.PLT_TYP:
        if   "pcolormesh" in my_fluid.PLT_TYP:
            X,Y = np.mgrid[ 1:my_fluid.NX-1, 1:my_fluid.NY-1 ]
            plt.pcolormesh(X,Y, the_data[3][1:-1, 1:-1], vmin=0,vmax=0.05)
            plt.colorbar()

            plt.quiver(X,Y, the_data[1][1:-1, 1:-1],
                            the_data[2][1:-1, 1:-1], color='w')

            plt_title  = "Velocity Arrows on Density Colormesh at Cycle"
            plt_title += " " + str(the_data[0])
            
        elif "contour" in my_fluid.PLT_TYP:
            plt.contourf(X,Y, the_data[4], alpha=0.5)
            plt.colorbar()

            plt.contour(X,Y, the_data[4])

            plt.quiver(X[::2, ::2],Y[::2, ::2], the_data[1][::2, ::2],
                                                the_data[2][::2, ::2])
            plt_title  = "Velocity Arrows on Pressure Contours at Cycle"
            plt_title += " " + str(the_data[0]) 
        else:
            plt.quiver(X[::3, ::3],Y[::3, ::3], the_data[1][::3, ::3],
                                                the_data[2][::3, ::3])
        plt.axis([X.min(),X.max(), Y.min(),Y.max()])
    else:
        ax = fig.gca(projection="3d")
        
        if   "velocity" in my_fluid.PLT_TYP:
            zlbl = "velocity"
            var  = the_data[1]
            ax.set_zlim(1, 2.5)
        elif "pressure" in my_fluid.PLT_TYP:
            zlbl = "pressure"
            var  = the_data[4]
        else:
            print("""Sorry, only non-quiver plots of velocity and
                  pressure are currently available. Please edit the
                  value of my_fluid.PLT_TYP accordingly.""")
            sys.exit()
        # end if/elif/else
                  
        if   "surface"   in my_fluid.PLT_TYP:
            ax.plot_surface(X,Y, var, rstride=1,cstride=1,
                            cmap=cm.coolwarm, linewidth=0,
                            antialiased=False)
        elif "wireframe" in my_fluid.PLT_TYP:
            ax.plot_wireframe(X,Y, var)
        else:
            print("""Sorry, the only non-quiver plots currently
                  available are surface and wireframe. Please edit the
                  value of my_fluid.PLT_TYP accordingly.""")
            sys.exit()
        # end if/elif/else
        ax.set_xlim(my_fluid.XMIN, my_fluid.XMAX)
        ax.set_ylim(my_fluid.YMIN, my_fluid.YMAX)
        ax.set_zlabel(zlbl)
        ax.view_init(30, 225)
    # end if/else
    plt.xlabel("x"); plt.ylabel("y")
    plt.title(plt_title)
    plt.savefig(plt_name)

    plt.clf()
    os.chdir(cur_path + "/StateFiles")
# end for
        
## end the_plotter_2d.py
