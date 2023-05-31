## implicit_vortex_2d.py
## by Ryan Farber 23 January  2015
## Last modified: 07 March    2015
"""
The purpose of this program is to apply the_argo using Ash's
Fluid for Dummies advect method (tweaked) to advect a vortex.

Note: two initial Runtime Warning's occur because of singular
nature of vortex initialization and are normal.
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
    u   = np.zeros((my_fluid.NX,my_fluid.NY)) # x-component of velocity
    v   = np.zeros((my_fluid.NX,my_fluid.NY)) # y-component of velocity
    cn  = np.zeros((my_fluid.NX,my_fluid.NY)) # concentration (tracers)
    p   = np.zeros((my_fluid.NX,my_fluid.NY)) # pressure
    src = np.zeros((my_fluid.NX,my_fluid.NY)) # src term for poisson eqn

        ##Initialize vortex
    x0 = 10.0;  y0 = 10.0          # center of vortex
    r2 = (XX-x0)**2 + (YY-y0)**2   # r2 is radius

    u[1:-1,1:-1] = np.where( r2 == 0, 0,  (YY-y0)/r2 )
    v[1:-1,1:-1] = np.where( r2 == 0, 0, -(XX-x0)/r2 )

        ##place some amount of tracer density at center 4 zones
    cn[ (len(cn)-1)/2 - 1, (len(cn)-1)/2 - 1 ] = 0.1
    cn[ (len(cn)-1)/2    , (len(cn)-1)/2 - 1 ] = 0.1
    cn[ (len(cn)-1)/2 - 1, (len(cn)-1)/2     ] = 0.1
    cn[ (len(cn)-1)/2    , (len(cn)-1)/2     ] = 0.1

    ##Save the state of the initial condition of the fluid
    cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                                   my_fluid.LABEL, ext)
    save_state([cycles, u,v, cn, p,src, "NA","NA", "NA","NA"],
               file_name)
else:
    data_file = glob.glob("*" + str(my_fluid.cycle_start) + "*")
    if data_file == []:
        print("Error! my_fluid.cycle_start data file not found.")
        sys.exit()
    # end if
    data_file = data_file[0]
    the_data = cPickle.load(open(data_file))
    u = the_data[1]; v = the_data[2]; cn = the_data[3]
    p = the_data[4]; src = the_data[5]
    the_data = 0 # to save memory
# end if


##Solve!
for cycles in xrange(my_fluid.cycle_start, my_fluid.NT+1):
    u_old = u.copy()
    
    u[1:-1,1:-1] = my_fluid.nonlinear_advect_implicit_2d(u, u,v, XX,YY)
    v[1:-1,1:-1] = my_fluid.nonlinear_advect_implicit_2d(v, u_old,v, XX,
                                                                     YY)
    u[1:-1,1:-1] = my_fluid.diffuse_implicit_2d(u)
    v[1:-1,1:-1] = my_fluid.diffuse_implicit_2d(v)
    
    u  = my_fluid.set_bnds_fixed_2dXI(u, -u[  1,  : ])
    u  = my_fluid.set_bnds_fixed_2dXF(u, -u[ -2,  : ])
    u  = my_fluid.set_bnds_fixed_2dYI(u,  u[  :,  1 ])
    u  = my_fluid.set_bnds_fixed_2dYF(u,  u[  :, -2 ])
    
    v  = my_fluid.set_bnds_fixed_2dXI(v,  v[  1,  : ])
    v  = my_fluid.set_bnds_fixed_2dXF(v,  v[ -2,  : ])
    v  = my_fluid.set_bnds_fixed_2dYI(v, -v[  :,  1 ])
    v  = my_fluid.set_bnds_fixed_2dYF(v, -v[  :, -2 ])
    
    u  = my_fluid.set_bnds_zero_flux_corners_2d( u)
    v  = my_fluid.set_bnds_zero_flux_corners_2d( v)
    
    src[1:-1,1:-1] = my_fluid.calc_source_2d(u,v)

    for dummy_var in xrange(my_fluid.NI):
        p = my_fluid.relax_pressure_poisson_2d(p, src)
        
        p = my_fluid.set_bnds_fixed_2dXI(p,  p[  1,  : ])
        p = my_fluid.set_bnds_fixed_2dXF(p,  p[ -2,  : ])
        p = my_fluid.set_bnds_fixed_2dYI(p,  p[  :,  1 ])
        p = my_fluid.set_bnds_fixed_2dYF(p,  p[  :, -2 ])
        p = my_fluid.set_bnds_zero_flux_corners_2d(p)
    # end for

    u[1:-1,1:-1] -= my_fluid.apply_pressure_2dX(p, 1.0/my_fluid.RHO)
    v[1:-1,1:-1] -= my_fluid.apply_pressure_2dY(p, 1.0/my_fluid.RHO)

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

    cn[1:-1,1:-1] = my_fluid.nonlinear_advect_implicit_2d(cn, u,v,XX,YY)
    cn[1:-1,1:-1] = my_fluid.diffuse_implicit_2d(cn)

    cn = my_fluid.set_bnds_fixed_2dXI(cn, cn[  1,  : ])
    cn = my_fluid.set_bnds_fixed_2dXF(cn, cn[ -2,  : ])
    cn = my_fluid.set_bnds_fixed_2dYI(cn, cn[  :,  1 ])
    cn = my_fluid.set_bnds_fixed_2dYF(cn, cn[  :, -2 ])

    cn = my_fluid.set_bnds_zero_flux_corners_2d(cn)

    if (cycles in my_fluid.SAVE_SP) or (cycles % my_fluid.SAVE_FREQ == 0):
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                           my_fluid.LABEL, ext)
        save_state([cycles, u,v, cn, p,src, "NA","NA", "NA","NA"],
                   file_name)
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state([cycles, u,v, cn, p,src, "NA","NA", "NA","NA"],
           file_name )

## end implicit_vortex_2d.py
