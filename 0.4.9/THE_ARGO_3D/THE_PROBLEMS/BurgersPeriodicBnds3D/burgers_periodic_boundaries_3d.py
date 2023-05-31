## burgers_periodic_boundaries_3d.py
## by Ryan Farber 30 December 2014
## Last modified: 24 January  2015
"""
The purpose of this program is to apply the_argo to propagate a hat
function by burgers equation (nonlinear convection and diffusion)
in two dimensions with periodic boundaries.
"""
import numpy as np
import sys; sys.path.insert(0, "../../");import the_bnds_setter_3d as bs 
sys.path.insert(0, "../../..")    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))
NX = my_fluid.NX; NY = my_fluid.NY; NZ = my_fluid.NZ; DT = my_fluid.DT
DX = my_fluid.DX; DY = my_fluid.DY; DZ = my_fluid.DZ; NU = my_fluid.NU

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
    u_old = u.copy(); v_old = v.copy(); w_old = w.copy()
    old_vars = [u_old, v_old, w_old]

    u = my_fluid.nonlinear_convect_3dX(u, u_old, v_old, w_old)
    v = my_fluid.nonlinear_convect_3dY(v, u_old, v_old, w_old)
    w = my_fluid.nonlinear_convect_3dZ(w, u_old, v_old, w_old)

    i = -1
    for var in [u,v,w]:
        i += 1
        var = my_fluid.diffuse_3d(var, old_vars[i])
        
        var = bs.set_bnds_periodic_diffusion_3dX(DT,DX,DY,DZ,NU,
                                                 var,old_vars[i])
        var = bs.set_bnds_periodic_diffusion_3dY(DT,DX,DY,DZ,NU,
                                                 var,old_vars[i])
        var = bs.set_bnds_periodic_diffusion_3dZ(DT,DX,DY,DZ,NU,
                                                 var,old_vars[i])

        var = bs.set_bnds_periodic_3dX(var)
        var = bs.set_bnds_periodic_3dY(var)
        var = bs.set_bnds_periodic_3dZ(var)

    if (cycles % my_fluid.SAVE_FREQ == 0):
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
        save_state( [cycles, u,v,w, 'NA','NA','NA','NA'], file_name ) 
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, u,v,w, 'NA','NA','NA','NA'], file_name )


## end burgers_periodic_boundaries_3d.py
