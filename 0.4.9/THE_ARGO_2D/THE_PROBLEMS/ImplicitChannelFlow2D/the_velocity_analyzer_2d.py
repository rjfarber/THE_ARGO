## the_velocity_analyzer_2d.py
## by Ryan Farber 17 February 2015
## Last modified: 28 February 2015
"""
The purpose of this program is to load saved variables calculated by the
Argo 2D and calculate the divergence of the velocity and calculate the
maximum of the absolute value of the flow variables.

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
import sys; sys.path.insert(0,"../../..")
from the_fluid            import The_Fluid
import os; cur_path = os.getcwd(); os.chdir(cur_path + "/StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Parameters
data_file_ext = ".p"
template = "{0:15}{1:15}{2:15}{3:15}{4:15}\n"

os.chdir(cur_path)
fid_write = open("velocity_analysis.txt", "w")
os.chdir(cur_path + "/StateFiles")

pfiles = glob.glob("*cycle*.p")
fid_write.write(template.format("Cycles","abs max u","abs max v",
                                "abs max cn","divV"))
for pfile in pfiles:
    the_data = cPickle.load(open(pfile, "rb"))

    divV  = (the_data[1][ 2:  , 1:-1 ] - the_data[1][  :-2, 1:-1 ]) \
                                      / (2*my_fluid.DX)
    divV += (the_data[2][ 1:-1, 2:   ] - the_data[2][ 1:-1,  :-2 ]) \
                                      / (2*my_fluid.DY)
    divV = np.absolute(divV)
    divV = sum(sum(divV))

    fid_write.write(template.format(str(the_data[0]),
str(np.max(np.absolute(the_data[1])))[:10],
str(np.max(np.absolute(the_data[2])))[:10],
str(np.max(np.absolute(the_data[3])))[:10],
str(divV)[:10]))

    os.chdir(cur_path + "/StateFiles")
# end for


## the_velocity_analyzer_2d.py
