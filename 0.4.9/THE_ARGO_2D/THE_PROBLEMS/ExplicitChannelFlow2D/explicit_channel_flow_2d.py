## explicit_channel_flow_2d.py
## by Ryan Farber 30 December 2014
## Last modified: 07 March    2015
"""
The purpose of this program is to apply the_argo to solve the
navier stokes equations in two dimensions for channel flow.

NOTE: After the first cycle, all the work has been done since
at that point the velocity has reached a steady state equilibrium;
hence the pressure is constant and, since we chose pressure = 0
initially, the constant end result is zero so don't plot the pressure
since matplotlib doesn't like it when pressure is all zero.
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
u   = np.zeros( (my_fluid.NX,my_fluid.NY) )# x-component of velocity
v   = np.zeros( (my_fluid.NX,my_fluid.NY) )# y-component of velocity
F   = np.ones(  (my_fluid.NX,my_fluid.NY) )# applied force (in x only)
p   = np.zeros( (my_fluid.NX,my_fluid.NY) )# pressure
src = np.zeros( (my_fluid.NX,my_fluid.NY) )# source term for poisson eqn

##Solve!
for cycles in xrange(1, my_fluid.NT+1):
    pct_err_vel = 1.0; stepcount = 0

    while pct_err_vel > my_fluid.MAE:
        stepcount += 1; u_old = u.copy()

        u[1:-1,1:-1] += (-my_fluid.nonlinear_advect_explicit_2d(u, u,v)
                      +   my_fluid.diffuse_explicit_2d(u)
                      +   my_fluid.apply_force_2d(F) )
        
        v[1:-1,1:-1] += (-my_fluid.nonlinear_advect_explicit_2d(v, u_old,v)
                      +   my_fluid.diffuse_explicit_2d(v) )
        
        ##Update ghost zones so boundaries are periodic (in X only)
        u[ 0, : ] = u[ -2,  : ]; u[ -1,  : ] = u[ 1, : ]
        v[ 0, : ] = v[ -2,  : ]; v[ -1,  : ] = v[ 1, : ]
        
        ##fixed along y
        u = my_fluid.set_bnds_fixed_2dYI(u,0)
        u = my_fluid.set_bnds_fixed_2dYF(u,0)
        v = my_fluid.set_bnds_fixed_2dYI(v,0)
        v = my_fluid.set_bnds_fixed_2dYF(v,0)
        
        src[ 1:-1, 1:-1 ] = my_fluid.calc_source_2d(u,v)

        for dummy_var in xrange(my_fluid.NI):
            p = my_fluid.relax_pressure_poisson_2d(p, src)

            ##Update ghost zones so boundaries are periodic (in X only)
            p[ 0, : ] = p[ -2,  : ]; p[ -1,  : ] = p[ 1, : ]
            ##dp/dy = 0 @ y = 0,2
            p[ :, 0 ] = p[ :, 1 ]; p[ :, -1 ] = p[ :, -2 ]
        # end for

        u[1:-1,1:-1] -= my_fluid.apply_pressure_2dX(p, 1.0/my_fluid.RHO)
        v[1:-1,1:-1] -= my_fluid.apply_pressure_2dY(p, 1.0/my_fluid.RHO)
                
        ##Update ghost zones so boundaries are periodic (in X only)
        u[ 0, : ] = u[ -2,  : ]; u[ -1,  : ] = u[ 1, : ]
        v[ 0, : ] = v[ -2,  : ]; v[ -1,  : ] = v[ 1, : ]
        
        ##fixed along y
        u = my_fluid.set_bnds_fixed_2dYI(u,0)
        u = my_fluid.set_bnds_fixed_2dYF(u,0)
        v = my_fluid.set_bnds_fixed_2dYI(v,0)
        v = my_fluid.set_bnds_fixed_2dYF(v,0)
 
        pct_err_vel = (np.sum(u) - np.sum(u_old)) / np.sum(u)
    # end while
    print(stepcount)
    
    file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                       my_fluid.LABEL, ext)
    save_state([cycles, u,v, "NA", p,src, F,"NA", "NA","NA"],
               file_name )
# end for


## end explicit_channel_flow_2d.py
