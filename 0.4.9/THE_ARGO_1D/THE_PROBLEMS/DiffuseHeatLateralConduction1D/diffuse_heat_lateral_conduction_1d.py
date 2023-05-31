## diffuse_heat_lateral_conduction_1d.py
## by Ryan Farber 28 December 2014
## Last modified: 26 January  2015
"""
The purpose of this program is to apply the_argo to find the steady
state solution of the IVBP of a heat equation of a rod, including heat
loss laterally and with each end fixed at 1, for an initial heat
distribution of one in the interior zones.
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
X = np.linspace(my_fluid.XMIN, my_fluid.XMAX, my_fluid.NX)
T = np.ones(my_fluid.NX)


##Save the state of the initial condition of the fluid
cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)

save_state( [cycles, T], file_name )


##Solve!
for cycles in xrange(1, my_fluid.NT+1):
    T = my_fluid.diffuse_1d(T) - (my_fluid.BETA * T)
    T[0] = 1; T[-1] = 1

    if (cycles % my_fluid.SAVE_FREQ == 0):
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
        save_state( [cycles, T], file_name)
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)

##Save State!
save_state( [cycles, T], file_name )


## end diffuse_heat_lateral_conduction_1d.py
