## double_vortex_opposite_sign_2d.py
## by Ryan Farber 27 January 2015
## Last modified: 03 May     2015
"""
The purpose of this program is to advect a double vortex.

Note: two initial Runtime Warning's occur because of the
singular nature of vortex initialization and are normal.
"""
import numpy as np
import sys; sys.path.insert(0, "../"); sys.path.insert(0, "../../")
sys.path.insert(0, "../../..") 
from the_bnds_setter_2d import The_Bnds_Setter_2D    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Setup
ext = ".p"  # filename extension for pickling
u   = np.zeros(  (my_fluid.NX,my_fluid.NY) ) # x-component of velocity
v   = np.zeros(  (my_fluid.NX,my_fluid.NY) ) # y-component of velocity
cn = np.zeros(  (my_fluid.NX,my_fluid.NY) ) # tracer particle density
p   = np.zeros(  (my_fluid.NX,my_fluid.NY) ) # pressure
src = np.zeros(  (my_fluid.NX,my_fluid.NY) ) # src term for poisson eqn

X,Y = np.mgrid[ 1:my_fluid.NX-1, 1:my_fluid.NY-1 ] # imitates 2d for loop

    ##Initialize double vortex
x0 = 25.0;  y0 = 5.0         # center of vortex
r2 = (X-x0)**2 + (Y-y0)**2   # r2 is radius

u[1:-1,1:-1] = np.where( r2 == 0, 0,  (Y-y0)/r2 )
v[1:-1,1:-1] = np.where( r2 == 0, 0, -(X-x0)/r2 )

u   = my_fluid.set_bnds_fixed_2dXI(u, -u[  1,  : ])
u   = my_fluid.set_bnds_fixed_2dXF(u, -u[ -2,  : ])
u   = my_fluid.set_bnds_fixed_2dYI(u,  u[  :,  1 ])
u   = my_fluid.set_bnds_fixed_2dYF(u,  u[  :, -2 ])

v   = my_fluid.set_bnds_fixed_2dXI(v,  v[  1,  : ])
v   = my_fluid.set_bnds_fixed_2dXF(v,  v[ -2,  : ])
v   = my_fluid.set_bnds_fixed_2dYI(v, -v[  :,  1 ])
v   = my_fluid.set_bnds_fixed_2dYF(v, -v[  :, -2 ])

x0 = 25.0;  y0 = 45.0         # center of vortex
r2 = (X-x0)**2 + (Y-y0)**2   # r2 is radius

u[1:-1,1:-1] += np.where( r2 == 0, 0, -(Y-y0)/r2 )
v[1:-1,1:-1] += np.where( r2 == 0, 0,  (X-x0)/r2 )

u   = my_fluid.set_bnds_fixed_2dXI(u, -u[  1,  : ])
u   = my_fluid.set_bnds_fixed_2dXF(u, -u[ -2,  : ])
u   = my_fluid.set_bnds_fixed_2dYI(u,  u[  :,  1 ])
u   = my_fluid.set_bnds_fixed_2dYF(u,  u[  :, -2 ])

v   = my_fluid.set_bnds_fixed_2dXI(v,  v[  1,  : ])
v   = my_fluid.set_bnds_fixed_2dXF(v,  v[ -2,  : ])
v   = my_fluid.set_bnds_fixed_2dYI(v, -v[  :,  1 ])
v   = my_fluid.set_bnds_fixed_2dYF(v, -v[  :, -2 ])

    ##place some amount of tracer density at center 4 zones
cn[ (len(cn)-1)/2 - 1, (len(cn)-1)/2 - 1 ] = 0.1
cn[ (len(cn)-1)/2    , (len(cn)-1)/2 - 1 ] = 0.1
cn[ (len(cn)-1)/2 - 1, (len(cn)-1)/2     ] = 0.1
cn[ (len(cn)-1)/2    , (len(cn)-1)/2     ] = 0.1


##Save the state of the initial condition of the fluid
cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, u,v, cn, p,src,"NA","NA","NA","NA"], file_name )


##Solve!
for cycles in xrange(1, my_fluid.NT+1):
        ##prep for new iteration
    p_old = p.copy()

        ##advect
    u[1:-1,1:-1]   = my_fluid.nonlinear_advect_implicit_2d( u, u,v, X,Y)
    v[1:-1,1:-1]   = my_fluid.nonlinear_advect_implicit_2d( v, u,v, X,Y)
    cn[1:-1,1:-1]  = my_fluid.nonlinear_advect_implicit_2d(cn, u,v, X,Y)    

        ##set zero flux boundary conditions
    u   = my_fluid.set_bnds_fixed_2dXI(u, -u[  1,  : ])
    u   = my_fluid.set_bnds_fixed_2dXF(u, -u[ -2,  : ])
    u   = my_fluid.set_bnds_fixed_2dYI(u,  u[  :,  1 ])
    u   = my_fluid.set_bnds_fixed_2dYF(u,  u[  :, -2 ])
    
    v   = my_fluid.set_bnds_fixed_2dXI(v,  v[  1,  : ])
    v   = my_fluid.set_bnds_fixed_2dXF(v,  v[ -2,  : ])
    v   = my_fluid.set_bnds_fixed_2dYI(v, -v[  :,  1 ])
    v   = my_fluid.set_bnds_fixed_2dYF(v, -v[  :, -2 ])

    cn = my_fluid.set_bnds_fixed_2dXI(cn,  cn[  1,  : ])
    cn = my_fluid.set_bnds_fixed_2dXF(cn,  cn[ -2,  : ])
    cn = my_fluid.set_bnds_fixed_2dYI(cn,  cn[  :,  1 ])
    cn = my_fluid.set_bnds_fixed_2dYF(cn,  cn[  :, -2 ])
    
    u  = my_fluid.set_bnds_zero_flux_corners_2d(u)
    v  = my_fluid.set_bnds_zero_flux_corners_2d(v)
    cn = my_fluid.set_bnds_zero_flux_corners_2d(cn)
    
        ##Begin Project
    src[1:-1,1:-1] = my_fluid.calc_source_2d(u,v)
    
    src = my_fluid.set_bnds_fixed_2dXI(src,  src[  1,  : ])
    src = my_fluid.set_bnds_fixed_2dXF(src,  src[ -2,  : ])
    src = my_fluid.set_bnds_fixed_2dYI(src,  src[  :,  1 ])
    src = my_fluid.set_bnds_fixed_2dYF(src,  src[  :, -2 ])

    src = my_fluid.set_bnds_zero_flux_corners_2d(src)

    for dummy_var in xrange(my_fluid.NI):
        p = my_fluid.relax_pressure_poisson_2d(p, src)
        
        p = my_fluid.set_bnds_fixed_2dXI(p,  p[  1,  : ])
        p = my_fluid.set_bnds_fixed_2dXF(p,  p[ -2,  : ])
        p = my_fluid.set_bnds_fixed_2dYI(p,  p[  :,  1 ])
        p = my_fluid.set_bnds_fixed_2dYF(p,  p[  :, -2 ])
        p = my_fluid.set_bnds_zero_flux_corners_2d(p)
    # end for

    u[1:-1,1:-1] -= my_fluid.apply_pressure_2dX(p, my_fluid.RHO)
    v[1:-1,1:-1] -= my_fluid.apply_pressure_2dY(p, my_fluid.RHO)

    u   = my_fluid.set_bnds_fixed_2dXI(u, -u[  1,  : ])
    u   = my_fluid.set_bnds_fixed_2dXF(u, -u[ -2,  : ])
    u   = my_fluid.set_bnds_fixed_2dYI(u,  u[  :,  1 ])
    u   = my_fluid.set_bnds_fixed_2dYF(u,  u[  :, -2 ])
    
    v   = my_fluid.set_bnds_fixed_2dXI(v,  v[  1,  : ])
    v   = my_fluid.set_bnds_fixed_2dXF(v,  v[ -2,  : ])
    v   = my_fluid.set_bnds_fixed_2dYI(v, -v[  :,  1 ])
    v   = my_fluid.set_bnds_fixed_2dYF(v, -v[  :, -2 ])

    u = my_fluid.set_bnds_zero_flux_corners_2d(u)
    v = my_fluid.set_bnds_zero_flux_corners_2d(v)
        ##End Project

    if (cycles % my_fluid.SAVE_FREQ == 0):
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
        save_state( [cycles, u,v, cn, p,src,'NA'], file_name )
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, u,v, cn, p,src,'NA'], file_name )


## end double_vortex_opposite_sign_2d.py
