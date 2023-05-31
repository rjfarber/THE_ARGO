## diffusion_periodic_boundaries_3d.py
## by Ryan Farber  9 January 2015
## Last modified: 24 Januray 2015
"""
The purpose of this program is to apply the_argo to propagate a hat
function by diffusion in three dimensions with periodic boundaries.
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
the_vars = [u,v,w]
old_vars = [u.copy(), v.copy(), w.copy()]
for cycles in xrange(1, my_fluid.NT+1):
    for i in xrange(len(the_vars)):
        old_vars[i] = the_vars[i].copy()
    # end for
    for i in xrange(len(the_vars)):
        the_vars[i] = my_fluid.diffuse_3d(the_vars[i], old_vars[i])
        
        the_vars[i] = bs.set_bnds_periodic_diffusion_3dX(
                    DT,DX,DY,DZ,NU, the_vars[i], old_vars[i])
        the_vars[i] = bs.set_bnds_periodic_diffusion_3dY(
                    DT,DX,DY,DZ,NU, the_vars[i], old_vars[i])
        the_vars[i] = bs.set_bnds_periodic_diffusion_3dZ(
                    DT,DX,DY,DZ,NU, the_vars[i], old_vars[i])
        
        the_vars[i] = bs.set_bnds_periodic_3dX(the_vars[i])
        the_vars[i] = bs.set_bnds_periodic_3dY(the_vars[i])
        the_vars[i] = bs.set_bnds_periodic_3dZ(the_vars[i])
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


## end diffusion_periodic_boundaries_3d.py
