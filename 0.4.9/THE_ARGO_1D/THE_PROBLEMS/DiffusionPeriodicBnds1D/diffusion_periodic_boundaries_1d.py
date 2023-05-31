## diffusion_periodic_boundaries_1d.py
## by Ryan Farber 29 December 2014
## Last modified: 22 January 2015
"""
The purpose of this program is to apply the_argo to propagate a hat
function by diffusion in one dimension with periodic boundaries.
"""
import numpy as np
import sys; sys.path.insert(0, "../..")
import the_bnds_setter_1d as bs
sys.path.insert(0, "../../..")
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
    u = bs.set_bnds_periodic_diffusion_1d(my_fluid.DT, my_fluid.DX,
                                          my_fluid.NU, u)
    u = bs.set_bnds_fixed_1dI(u, u[-1])

    if (cycles % my_fluid.SAVE_FREQ == 0):
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)       
        save_state( [cycles, u], file_name )
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)

##Save State!
save_state( [cycles, u], file_name )


## end diffusion_periodic_boundaries_1d.py
