## linear_advection_periodic_1d.py
## by Ryan Farber 29 December 2014
## Last modified: 20 December 2015
"""
The purpose of this program is to apply the_argo to propagate a hat
function by linear advection in one dimension with periodic boundaries.
"""
import numpy as np
import sys; sys.path.insert(0, "../../"); sys.path.insert(0, "../../..")
from the_fluid import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))


##Setup
u = np.ones(my_fluid.NX)
u[ 0.25/my_fluid.DX : 0.5/my_fluid.DX+1 ] = 2.0


##Save the state of the initial condition of the fluid
cycles = 1; file_name = my_fluid.get_file_name(cycles,
            my_fluid.MAX_CYC_MAG, my_fluid.LABEL, my_fluid.OUT_FILE_EXT)
my_fluid.save_state( [cycles, u], file_name )


##Solve!
for cycles in xrange(1, my_fluid.NT+1):

    u = my_fluid.linear_advect_1d(u)
    
    u[0] = u[-2]; u[-1] = u[1]

    if (cycles % my_fluid.SAVE_FREQ == 0):
        file_name = my_fluid.get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                    my_fluid.LABEL, my_fluid.OUT_FILE_EXT)
        my_fluid.save_state( [cycles, u], file_name ) 
    # end if
# end for
file_name = my_fluid.get_file_name(cycles, my_fluid.MAX_CYC_MAG,
            my_fluid.LABEL, my_fluid.OUT_FILE_EXT)

##Save State!
my_fluid.save_state( [cycles, u], file_name )


## end linear_advection_periodic_1d.py
