## the_div_calcer_2d.py
## by Ryan Farber 23 January  2015
## Last modified: 15 February 2015
"""
The purpose of this program is to load velocities calculated by the
Argo 2D and calculate the divergence of the velocity field.

  NOTE: This file must be copied by a bash script to the specific
Problem file folder it shall be used in to work.

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


os.chdir(cur_path)
fid_write = open("divergence.txt", "w")
os.chdir(cur_path + "/StateFiles")

pfiles = glob.glob("*cycle*.p")
for pfile in pfiles:
    the_data = cPickle.load(open(pfile, "rb"))

    divV  = (the_data[1][ 2:  , 1:-1 ] - the_data[1][  :-2, 1:-1 ]) \
                                      / (2*my_fluid.DX)
    divV += (the_data[2][ 1:-1, 2:   ] - the_data[2][ 1:-1,  :-2 ]) \
                                      / (2*my_fluid.DY)
    divV = np.absolute(divV)
    divV = sum(sum(divV))

    fid_write.write("\nCycle: " + str(the_data[0]) +
                    "; divV: " + str(divV))
    if len(the_data) == 8:
        if the_data[-1] != "NA":
            divB  = (the_data[6][ 2:  , 1:-1 ]
                   - the_data[6][  :-2, 1:-1 ]) / (2*my_fluid.DX)
            divB += (the_data[7][ 1:-1, 2:   ]
                   - the_data[7][ 1:-1,  :-2 ]) / (2*my_fluid.DY)
            divB = np.absolute(divB)
            divB = sum(sum(divB))

            fid_write.write("\tdivB: " + str(divB))
        # end if
    # end if
    os.chdir(cur_path + "/StateFiles")
# end for


## the_div_calcer_2d.py
