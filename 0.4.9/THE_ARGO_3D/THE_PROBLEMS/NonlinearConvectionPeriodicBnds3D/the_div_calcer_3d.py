## the_div_calcer_3d.py
## by Ryan Farber 26 January 2015
## Last modified: 27 January 2015
"""
The purpose of this program is to load velocities calculated by the
Argo 3D and calculate the divergence of the velocity field.

  NOTE: This file must be copied by a bash script to the specific
Problem file folder it shall be used in to work.

Indexing of "the_data" depends on what dimension of "save_state"
was used; the below disambiguates for data from The Argo 3D.

the_data[0] = cycles
the_data[1] = u       # x-component of velocity
the_data[2] = v       # y-component of velocity
the_data[3] = w       # z-component of velocity
the_data[4] = rho     # tracer density
the_data[5] = p       # pressure
the_data[6] = src     # src term of pressure poisson eqn
the_data[7] = F       # external force
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

    div  = (the_data[1][ 2:  , 1:-1, 1:-1 ]
         -  the_data[1][  :-2, 1:-1, 1:-1 ]) / (2*my_fluid.DX)
    
    div += (the_data[2][ 1:-1, 2:  , 1:-1 ]
         -  the_data[2][ 1:-1,  :-2, 1:-1 ]) / (2*my_fluid.DY)

    div += (the_data[3][ 1:-1, 1:-1, 2:   ]
         -  the_data[3][ 1:-1, 1:-1,  :-2 ]) / (2*my_fluid.DZ)
    
    div = np.absolute(div)
    div = sum(sum(sum(div)))

    fid_write.write("Cycle: " + str(the_data[0]) +
                    "; Divergence: " + str(div) + "\n")
    os.chdir(cur_path + "/StateFiles")
# end for
        
## the_div_calcer_3d.py
