## channel_flow_div_only_3d.py
## by Ryan Farber 27 January 2015
## Last modified: 27 January 2015
"""
The purpose of this program is to apply the_argo to solve the
navier stokes equations in three dimensions for channel flow,
including only the divergence terms in calculating the source for the
pressure poisson equation.

NOTE: After the first cycle, all the work has been done since
at that point the velocity has reached a steady state equilibrium;
hence the pressure is constant and, since we chose pressure = 0
initially, the constant end result is zero so don't do
plot_full_3d since it doesn't like it when pressure is all zero.

Questionable if this really works correct.
"""
import numpy as np
import sys; sys.path.insert(0, "../../");import the_bnds_setter_3d as bs 
sys.path.insert(0, "../../..")    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))
NX = my_fluid.NX; NY  = my_fluid.NY; NZ = my_fluid.NZ; DT  = my_fluid.DT
DX = my_fluid.DX; DY  = my_fluid.DY; DZ = my_fluid.DZ
NU = my_fluid.NU; RHO = my_fluid.RHO

##Setup
ext = ".p"  # filename extension for pickling
u   = np.zeros( (NX,NY,NZ) ) # x-component of velocity
v   = np.zeros( (NX,NY,NZ) ) # y-component of velocity
w   = np.zeros( (NX,NY,NZ) ) # z-component of velocity
F   = np.ones(  (NX,NY,NZ) ) # applied force
p   = np.zeros( (NX,NY,NZ) ) # pressure
src = np.zeros( (NX,NY,NZ) ) # source term for poisson eqn


##Solve!
for cycles in xrange(1, my_fluid.NT+1):
    pct_err_vel = 1.0; stepcount = 0

    while pct_err_vel > my_fluid.MAE:
        stepcount += 1
        
            ##prep pressure
        src = my_fluid.calc_source_div_only_3d(src, u,v,w)
        src = bs.set_bnds_periodic_src_div_only_3dX(DT,DX,DY,DZ,RHO,src,
                                                   u,v,w)

            ##calc pressure
        for dummy_var in xrange(my_fluid.NI):
            p_old = p.copy()

            p = my_fluid.relax_pressure_poisson_3d(p,p_old, src)
            p = bs.set_bnds_periodic_pressure_3dX(DX,DY,DZ, p,p_old,src)

            ##Boundary Conditions
                ##dp/dy = 0 @ y = 0,2
            p[ :, 0, : ] = p[ :, 1, : ]; p[ :, -1, : ] = p[ :, -2, : ]
        # end for

            ##prep velocity
        u_old = u.copy()
        v_old = v.copy()
        w_old = w.copy()
        
            ##calc velocity
        u = my_fluid.nonlinear_convect_3dX(u, u_old, v_old, w_old)
        v = my_fluid.nonlinear_convect_3dY(v, u_old, v_old, w_old)
        w = my_fluid.nonlinear_convect_3dZ(w, u_old, v_old, w_old)
        
        u = my_fluid.diffuse_3d(u, u_old)
        v = my_fluid.diffuse_3d(v, v_old)
        w = my_fluid.diffuse_3d(w, w_old)
        
        u = my_fluid.apply_pressure_3dX(u, p)
        v = my_fluid.apply_pressure_3dY(v, p)
        w = my_fluid.apply_pressure_3dZ(w, p)
        
        u = my_fluid.apply_force_3d(u, F)

        ##boundary conditions
            ##periodic along x
        u = bs.set_bnds_periodic_diffusion_3dX(DT,DX,DY,DZ,NU,
                                                     u,u_old)
        v = bs.set_bnds_periodic_diffusion_3dX(DT,DX,DY,DZ,NU,
                                                     v,v_old)
        w = bs.set_bnds_periodic_diffusion_3dX(DT,DX,DY,DZ,NU,
                                                     w,w_old)
        
        u,v,w = bs.set_bnds_periodic_apply_pressure_3dX(DT,DX,DY,
                                                    DZ,RHO,u,v,w, p)
        u = bs.set_bnds_periodic_apply_force_3dX(DT, u, F)
        
        u = bs.set_bnds_periodic_3dX(u)
        v = bs.set_bnds_periodic_3dX(v)
        w = bs.set_bnds_periodic_3dX(w)
        
            ##fixed along y
        u[ :, 0, : ] = 0; u[ :, -1,  : ] = 0
        v[ :, 0, : ] = 0; v[ :, -1,  : ] = 0
        w[ :, 0, : ] = 0; w[ :, -1,  : ] = 0

            ##fixed along z
        u[ :, :, 0 ] = 0; u[ :,  :, -1 ] = 0
        v[ :, :, 0 ] = 0; v[ :,  :, -1 ] = 0
        w[ :, :, 0 ] = 0; w[ :,  :, -1 ] = 0

        pct_err_vel = ((np.sum(u) - np.sum(u_old))
                    /   np.sum(u))
    # end while
    print(stepcount)
    
    file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
    save_state( [cycles, u,v,w, 'NA' ,p,src,F], file_name )
# end for


## end channel_flow_div_only_3d.py
