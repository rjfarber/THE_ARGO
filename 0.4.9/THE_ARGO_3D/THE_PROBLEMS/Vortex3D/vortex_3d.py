## vortex_3d.py
## by Ryan Farber 16 January 2015
## Last modified: 24 January 2015
"""
The purpose of this program is to apply the_argo using Ash's
Fluid for Dummies advect method to advect a vortex.
    Note: two initial Runtime Warning's occur because of singular
nature of vortex initialization and are normal.
"""
import numpy as np
import sys; sys.path.insert(0, "../../");import the_bnds_setter_3d as bs 
sys.path.insert(0, "../../..")    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))
N = my_fluid.NX; DT = my_fluid.DT

##Setup
ext = ".p"  # filename extension for pickling
u       = np.zeros( (N,N,N) ) # x-component of velocity
v       = np.zeros( (N,N,N) ) # y-component of velocity
w       = np.zeros( (N,N,N) ) # z-component of velocity
rho     = np.zeros( (N,N,N) ) # tracer particle density
p       = np.zeros( (N,N,N) ) # pressure
src     = np.zeros( (N,N,N) ) # src term for poisson pressure eqn

X,Y,Z = np.mgrid[ 1:N-1, 1:N-1, 1:N-1 ] # imitates triple for loop

    ##Initialize vortex
x0 = 10.0;  y0 = 5.0         # center of vortex
r2 = (X-x0)**2 + (Y-y0)**2   # r2 is radius

u[ 1:-1, 1:-1, 1:-1 ] = np.where( r2 == 0, 0,  (Y-y0)/r2 )
v[ 1:-1, 1:-1, 1:-1 ] = np.where( r2 == 0, 0, -(X-x0)/r2 )

    ##place some amount of tracer density at center 8 zones
rho[ N/2-1, N/2-1, N/2-1 ] = 0.1
rho[ N/2  , N/2-1, N/2-1 ] = 0.1
rho[ N/2-1, N/2  , N/2-1 ] = 0.1
rho[ N/2-1, N/2-1, N/2   ] = 0.1
rho[ N/2  , N/2  , N/2-1 ] = 0.1
rho[ N/2  , N/2-1, N/2   ] = 0.1
rho[ N/2-1, N/2  , N/2   ] = 0.1
rho[ N/2  , N/2  , N/2   ] = 0.1


##Save the state of the initial condition of the fluid
cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, u,v,w, rho,p,src, 'NA'], file_name )


##Solve!
for cycles in xrange(1, my_fluid.NT+1):
        ##prep for new iteration
    p_old = p.copy()

        ##prep advect (velocity and density)
    x = X - (DT * u[ 1:-1, 1:-1, 1:-1 ])
    y = Y - (DT * v[ 1:-1, 1:-1, 1:-1 ])
    z = Z - (DT * w[ 1:-1, 1:-1, 1:-1 ])

    x = np.where(x < 0.5, 0.5, x)
    y = np.where(y < 0.5, 0.5, y)
    z = np.where(z < 0.5, 0.5, z)
    
    x = np.where(x > (N-2) + 0.5, (N-2) + 0.5, x)
    y = np.where(y > (N-2) + 0.5, (N-2) + 0.5, y)
    z = np.where(z > (N-2) + 0.5, (N-2) + 0.5, z)

    i0 = x.astype(int); j0 = y.astype(int); k0 = z.astype(int)
    i1 = i0 + 1; j1 = j0 + 1; k1 = k0 + 1

    s1 = x - i0; t1 = y - j0; u1 = z - k0
    s0 = 1 - s1; t0 = 1 - t1; u0 = 1 - u1

        ##advect
    u   =my_fluid.advect_3d(u, i0,i1, j0,j1, k0,k1, s0,s1, t0,t1, u0,u1)
    v   =my_fluid.advect_3d(v, i0,i1, j0,j1, k0,k1, s0,s1, t0,t1, u0,u1)
    w   =my_fluid.advect_3d(w, i0,i1, j0,j1, k0,k1, s0,s1, t0,t1, u0,u1)
    rho=my_fluid.advect_3d(rho,i0,i1, j0,j1, k0,k1, s0,s1, t0,t1, u0,u1)    

        ##set zero flux boundary conditions
    u   = bs.set_bnds_zero_flux_3dX(u)
    v   = bs.set_bnds_zero_flux_3dY(v)
    w   = bs.set_bnds_zero_flux_3dZ(w)
    rho = bs.set_bnds_zero_flux_scalar_3d(rho)

    u   = bs.set_bnds_zero_flux_corners_3d(u)
    v   = bs.set_bnds_zero_flux_corners_3d(v)
    w   = bs.set_bnds_zero_flux_corners_3d(w)
    rho = bs.set_bnds_zero_flux_corners_3d(rho)
    
        ##Begin Project
    src = my_fluid.calc_source_3d(src, u,v,w)

    src = bs.set_bnds_zero_flux_scalar_3d(src)
    src = bs.set_bnds_zero_flux_corners_3d(src)

    for dummy_var in xrange(my_fluid.NI):
        p = my_fluid.relax_pressure_poisson_3d(p,p_old, src)
        p = bs.set_bnds_zero_flux_scalar_3d(p)
        p = bs.set_bnds_zero_flux_corners_3d(p)
    # end for

    u = my_fluid.apply_pressure_3dX(u, p)
    v = my_fluid.apply_pressure_3dY(v, p)
    w = my_fluid.apply_pressure_3dZ(w, p)

    u = bs.set_bnds_zero_flux_3dX(u)
    v = bs.set_bnds_zero_flux_3dY(v)
    w = bs.set_bnds_zero_flux_3dZ(w)

    u = bs.set_bnds_zero_flux_corners_3d(u)
    v = bs.set_bnds_zero_flux_corners_3d(v)
    w = bs.set_bnds_zero_flux_corners_3d(w)
        ##End Project

    if (cycles % my_fluid.SAVE_FREQ == 0):
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
        save_state( [cycles, u,v,w, rho,p,src,'NA'], file_name ) 
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, u,v,w, rho,p,src,'NA'], file_name )


## end vortex_3d.py
