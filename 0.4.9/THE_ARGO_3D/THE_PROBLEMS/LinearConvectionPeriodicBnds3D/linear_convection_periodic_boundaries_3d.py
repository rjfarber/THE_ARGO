## linear_convecton_periodic_boundaries_3d.py
## by Ryan Farber  8 January 2015
## Last modified: 24 January 2015
"""
The purpose of this program is to apply the_argo to propagate a hat
function by linear convection in three dimensions
with periodic boundaries.
"""
import numpy as np
import sys; sys.path.insert(0, "../../");import the_bnds_setter_3d as bs 
sys.path.insert(0, "../../..")    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))
NX = my_fluid.NX; NY = my_fluid.NY; NZ = my_fluid.NZ
DX = my_fluid.DX; DY = my_fluid.DY

##Setup
ext = ".p"  # filename extension for pickling
u = np.ones( (NX,NY,NZ) ) # x-component of velocity
v = np.ones( (NX,NY,NZ) ) # y-component of velocity
w = np.ones( (NX,NY,NZ) ) # z-component of velocity

u[ 0.5/DX : 1.0/DX+1, 0.5/DY : 1.0/DY+1, : ] = 2.0
v[ 0.5/DX : 1.0/DX+1, 0.5/DY : 1.0/DY+1, : ] = 2.0
w[ 0.5/DX : 1.0/DX+1, 0.5/DY : 1.0/DY+1, : ] = 2.0


##Save the state of the initial condition of the fluid
cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, u,v,w, 'NA','NA','NA','NA'], file_name )


##Solve!
for cycles in xrange(1, my_fluid.NT+1):
    u_old = u.copy()
    v_old = v.copy()
    w_old = w.copy()
    
    u = my_fluid.linear_convect_3d(u, u_old)
    v = my_fluid.linear_convect_3d(v, v_old)
    w = my_fluid.linear_convect_3d(w, w_old)

    for var in [u,v,w]:
        var = bs.set_bnds_periodic_3dX(var)
        var = bs.set_bnds_periodic_3dY(var)
        var = bs.set_bnds_periodic_3dZ(var)
    # end for
    
    if (cycles % my_fluid.SAVE_FREQ == 0):
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
        save_state( [cycles, u,v,w, 'NA','NA','NA','NA'], file_name ) 
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, u,v,w, 'NA','NA','NA','NA'], file_name )


## end linear_convection_periodic_boundaries_3d.py
