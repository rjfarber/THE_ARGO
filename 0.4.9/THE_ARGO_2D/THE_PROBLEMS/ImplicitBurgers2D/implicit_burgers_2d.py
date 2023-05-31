## implicit_burgers_2d.py
## by Ryan Farber 12 February 2015
## Last modified: 07 March    2015
"""
The purpose of this program is to apply the_argo to propagate a hat
function by burgers equation (nonlinear advection and diffusion)
in two dimensions.
"""
import glob
import numpy as np
import sys; sys.path.insert(0, "../../..")    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Setup
ext = ".p"  # filename extension for pickling
XX,YY = np.mgrid[1:my_fluid.NX-1, 1:my_fluid.NY-1]# imitates 2d for loop

if my_fluid.cycle_start == 1:
    u  = np.ones( (my_fluid.NX,my_fluid.NY) ) # x-component of velocity
    v  = np.ones( (my_fluid.NX,my_fluid.NY) ) # y-component of velocity

    u[ 0.5/my_fluid.DX : 1.0/my_fluid.DX+1,
       0.5/my_fluid.DY : 1.0/my_fluid.DY+1 ] = 2.0
    v[ 0.5/my_fluid.DX : 1.0/my_fluid.DX+1,
       0.5/my_fluid.DY : 1.0/my_fluid.DY+1 ] = 2.0

    ##Save the state of the initial condition of the fluid
    cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                                   my_fluid.LABEL, ext)
    save_state([cycles, u,v, "NA", "NA","NA", "NA","NA", "NA","NA"],
               file_name )
else:
    data_file = glob.glob("*" + str(my_fluid.cycle_start) + "*")
    if data_file == []:
        print("Error! my_fluid.cycle_start data file not found.")
        sys.exit()
    # end if
    data_file = data_file[0]
    the_data = cPickle.load(open(data_file))
    u = the_data[1]; v = the_data[2]
    the_data = 0 # to save memory
# end if


##Solve!
for cycles in xrange(my_fluid.cycle_start, my_fluid.NT+1):

    u[1:-1,1:-1] = my_fluid.nonlinear_advect_implicit_2d(u, u,v, XX,YY)
    v[1:-1,1:-1] = my_fluid.nonlinear_advect_implicit_2d(v, u,v, XX,YY)
    u[1:-1,1:-1] = my_fluid.diffuse_implicit_2d(u)
    v[1:-1,1:-1] = my_fluid.diffuse_implicit_2d(v)

    u = my_fluid.set_bnds_fixed_2d(u, 1)
    v = my_fluid.set_bnds_fixed_2d(v, 1)
    
    if (cycles % my_fluid.SAVE_FREQ == 0):
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                           my_fluid.LABEL, ext)
        save_state([cycles, u,v, "NA", "NA","NA", "NA","NA", "NA","NA"],
               file_name )
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                   my_fluid.LABEL, ext)
save_state([cycles, u,v, "NA", "NA","NA", "NA","NA", "NA","NA"],
               file_name )


## end implicit_burgers_2d.py
