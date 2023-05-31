## cavity_flow_3d.py
## by Ryan Farber 15 January 2015
## Last modified: 24 January 2015
"""
The purpose of this program is to apply the_argo to solve the
navier stokes equations in three dimensions for cavity flow.
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
u = np.zeros( (NX,NY,NZ) ) # x-component of velocity vector
v = np.zeros( (NX,NY,NZ) ) # y-component of velocity vector
w = np.zeros( (NX,NY,NZ) ) # z-component of velocity vector

    ##Boundary condition, fixed at one face (the 'cavity' or 'lid')
u[:,-1,:] = 1.0

p   = np.zeros( (NX,NY,NZ) ) # pressure
src = np.zeros( (NX,NY,NZ) ) # source term for poisson eqn


##Solve!
for cycles in xrange(1, my_fluid.NT+1):
        ##first pressure
    src = my_fluid.calc_source_3d(src, u,v,w)
    pct_err = 1.0; counter = 0

    while pct_err > my_fluid.MAE and counter < my_fluid.NI:
        ##prep for new iteration
        p_old = p.copy()
        counter += 1
        
        p = my_fluid.relax_pressure_poisson_3d(p,p_old, src)
        p = bs.set_bnds_pressure_3d(p, 'cavity_flow')

        ##calc percent error for this iteration
        pct_err = (np.sum(np.abs(p) - np.abs(p_old))
                /  np.sum(np.abs(p_old)))
    # end while
        ##now velocity
    u_old = u.copy()
    v_old = v.copy()
    w_old = w.copy()
    
    u = my_fluid.nonlinear_convect_3dX(u, u_old, v_old, w_old)
    v = my_fluid.nonlinear_convect_3dY(v, u_old, v_old, w_old)
    w = my_fluid.nonlinear_convect_3dZ(w, u_old, v_old, w_old)
    
    u = my_fluid.diffuse_3d(u, u_old)
    v = my_fluid.diffuse_3d(v, v_old)
    w = my_fluid.diffuse_3d(w, w_old)

    u = my_fluid.apply_pressure_3dX(u, p)
    v = my_fluid.apply_pressure_3dY(v, p)
    w = my_fluid.apply_pressure_3dZ(w, p)

    u = bs.set_bnds_fixed_3dX(u, 0, 0)
    u = bs.set_bnds_fixed_3dY(u, 0, 1)
    u = bs.set_bnds_fixed_3dZ(u, 0, 0)

    v = bs.set_bnds_fixed_3dX(v, 0, 0)
    v = bs.set_bnds_fixed_3dY(v, 0, 0)
    v = bs.set_bnds_fixed_3dZ(v, 0, 0)

    w = bs.set_bnds_fixed_3dX(w, 0, 0)
    w = bs.set_bnds_fixed_3dY(w, 0, 0)
    w = bs.set_bnds_fixed_3dZ(w, 0, 0)
    
    if cycles in [2, 200]:
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
        save_state( [cycles, u,v,w, 'NA',p,src,'NA'], file_name ) 
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, u,v,w, 'NA', p,src, 'NA'], file_name )


## end cavity_flow_3d.py
