## the_analyzer_2d.py
## by Ryan Farber 17 February 2015
## Last modified: 07 March    2015
"""
The purpose of this program is to load saved variables calculated by the
Argo 2D and calculate the divergence of the velocity and magnetic fields
and calculate the maximum of the absolute value of the flow variables.

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
import matplotlib.pyplot as plt
import glob
import numpy as np
import sys; sys.path.insert(0,"../../..")
from the_fluid            import The_Fluid
import os; cur_path = os.getcwd(); os.chdir(cur_path + "/StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Parameters
data_file_ext = ".p"
template = "{0:10}{1:15}{2:10}\t{3:10}\t{4:10}\t{5:10}\t  {6:10}\t  {7:10}\n"

os.chdir(cur_path)
fid_write = open("analysis.txt", "w")
os.chdir(cur_path + "/StateFiles")

pfiles = glob.glob("*cycle*.p")
fid_write.write(template.format("Cycles","abs max u","abs max v",
                                "abs max cn","abs max Bx","abs max By",
                                "divV", "divB"))
template = "{0}\t{1:10.2e}\t{2:10.2e}\t{3:10.2e}\t{4:10.2e}\t{5:10.2e}\t{6:10.2e}\t{7:10.2e}\n"

X = np.linspace(my_fluid.XMIN,my_fluid.XMAX,my_fluid.NX-2)
Y = np.linspace(my_fluid.YMIN,my_fluid.YMAX,my_fluid.NY-2)
Y,X = np.meshgrid(X,Y)

for pfile in pfiles:
    the_data = cPickle.load(open(pfile, "rb"))

    os.chdir(cur_path + "/Plots")

    divV  = (the_data[1][ 2:  , 1:-1 ] - the_data[1][  :-2, 1:-1 ]) \
                                      / (2*my_fluid.DX)
    divV += (the_data[2][ 1:-1, 2:   ] - the_data[2][ 1:-1,  :-2 ]) \
                                      / (2*my_fluid.DY)
##    plt.pcolormesh(X,Y, divV)
##    plt.title("divV")
##    plt.xlabel("x")
##    plt.ylabel("y")
##    plt.colorbar()
##    plt.savefig("divV_cycle_" + str(the_data[0]) + ".png")
##    plt.close()

    divV = np.absolute(divV)
    divV = sum(sum(divV))

    if len(the_data) >= 10:
        if the_data[-1] != "NA":
            divB  = (the_data[8][ 2:  , 1:-1 ]
                   - the_data[8][  :-2, 1:-1 ]) / (2*my_fluid.DX)
            divB += (the_data[9][ 1:-1, 2:   ]
                   - the_data[9][ 1:-1,  :-2 ]) / (2*my_fluid.DY)

##            plt.pcolormesh(X,Y, divB)
##            plt.title("divB")
##            plt.xlabel("x")
##            plt.ylabel("y")
##            plt.colorbar()
##            plt.savefig("divB_cycle_" + str(the_data[0]) + ".png")
##            plt.close()
            
            divB = np.absolute(divB)
            divB = sum(sum(divB))
        else:
            divB = "NA"
        # end if/else
    # end if
    fid_write.write(template.format(the_data[0],
np.max(np.absolute(the_data[1])),
np.max(np.absolute(the_data[2])),
np.max(np.absolute(the_data[3])),
np.max(np.absolute(the_data[8])),
np.max(np.absolute(the_data[9])), divV, divB))

    os.chdir(cur_path + "/StateFiles")
# end for


## the_analyzer_2d.py
