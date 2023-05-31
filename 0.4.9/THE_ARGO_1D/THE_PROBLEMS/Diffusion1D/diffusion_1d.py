## diffusion_1d.py
## by Ryan Farber 28 December 2014
## Last modified: 22 January 2015
"""
The purpose of this program is to apply the_argo to propagate a hat
function by diffusion in one dimension.
"""
import numpy as np
import sys; sys.path.insert(0, "../../..")
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))


##Setup
ext = ".p"  # filename extension for pickling
u = np.ones(my_fluid.NX)
u[ 0.5/my_fluid.DX : 1.0/my_fluid.DX+1 ] = 2.0


##Save the state of the initial condition of the fluid
cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)

save_state( [cycles, u], file_name )


##Solve!
for cycles in xrange(1, my_fluid.NT+1):
    u = my_fluid.diffuse_1d(u)
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)

##Save State!
save_state( [cycles, u], file_name )


## end diffusion_1d.py
