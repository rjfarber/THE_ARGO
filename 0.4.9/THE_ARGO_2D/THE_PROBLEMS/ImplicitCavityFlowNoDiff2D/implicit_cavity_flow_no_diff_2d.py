## implicit_cavity_flow_fft_2d.py
## by Ryan Farber 03 May 2015
## Last modified: 03 May 2015
"""
The purpose of this program is to apply the_argo to solve the
navier stokes equations in two dimensions for cavity flow.
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
    u = np.zeros( (my_fluid.NX,my_fluid.NY) ) # x-component of velocity
    v = np.zeros( (my_fluid.NX,my_fluid.NY) ) # y-component of velocity

        ##Boundary condition, fixed at one wall (the 'cavity' or 'lid')
    u[:,-1] = 1.0

    ##Pressure
    p   = np.zeros( (my_fluid.NX,my_fluid.NY) )# pressure
    src = np.zeros( (my_fluid.NX,my_fluid.NY) )# source term for poisson eqn
else:
    data_file = glob.glob("*" + str(my_fluid.cycle_start) + "*")
    if data_file == []:
        print("Error! my_fluid.cycle_start data file not found.")
        sys.exit()
    # end if
    data_file = data_file[0]
    the_data = cPickle.load(open(data_file))
    u = the_data[1]; v = the_data[2]; p = the_data[4]; src = the_data[5]
    the_data = 0 # to save memory
# end if


##Solve!
for cycles in xrange(my_fluid.cycle_start, my_fluid.NT+1):

    u[1:-1,1:-1] = my_fluid.nonlinear_advect_implicit_2d(u, u,v, XX,YY)
    v[1:-1,1:-1] = my_fluid.nonlinear_advect_implicit_2d(v, u,v, XX,YY)

    u = my_fluid.set_bnds_fixed_2d(u, 0)
    u = my_fluid.set_bnds_fixed_2dYF(u, 1)
    v = my_fluid.set_bnds_fixed_2d(v, 0)
    
    u[1:-1,1:-1] = my_fluid.diffuse_implicit_2d(u_old,u,my_fluid.NU)
    v[1:-1,1:-1] = my_fluid.diffuse_implicit_2d(v_old,v,my_fluid.NU)

    u = my_fluid.set_bnds_fixed_2d(u, 0)
    u = my_fluid.set_bnds_fixed_2dYF(u, 1)
    v = my_fluid.set_bnds_fixed_2d(v, 0) 
    
    divV  = (u[ 2:  , 1:-1 ] - u[  :-2, 1:-1 ]) \
                                      / (2*my_fluid.DX)
    divV += (v[ 1:-1, 2:   ] - v[ 1:-1,  :-2 ]) \
                                      / (2*my_fluid.DY)
    divV = sum(sum(np.absolute(divV)))

    while divV > 1e-3:
        src[1:-1,1:-1] = my_fluid.calc_source_2d(u,v)
                
        p = my_fluid.relax_pressure_poisson_2d(p, src)

        u[1:-1,1:-1] -= my_fluid.apply_pressure_2dX(p, 1.0/my_fluid.RHO)
        v[1:-1,1:-1] -= my_fluid.apply_pressure_2dY(p, 1.0/my_fluid.RHO)
        
        ##Update ghost zones
        u = my_fluid.set_bnds_fixed_2d(u, 0)
        u = my_fluid.set_bnds_fixed_2dYF(u, 1)
        v = my_fluid.set_bnds_fixed_2d(v, 0)

        divV  = (u[ 2:  , 1:-1 ] - u[  :-2, 1:-1 ]) \
                                      / (2*my_fluid.DX)
        divV += (v[ 1:-1, 2:   ] - v[ 1:-1,  :-2 ]) \
                                      / (2*my_fluid.DY)
        divV = sum(sum(np.absolute(divV)))
    # end while  

    if cycles in [2, 200]:
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                           my_fluid.LABEL, ext)
        save_state([cycles, u,v, "NA", p,src, "NA","NA", "NA","NA"],
                   file_name ) 
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                   my_fluid.LABEL, ext)
save_state([cycles, u,v, "NA", p,src, "NA","NA", "NA","NA"],
                   file_name )  


## end implicit_cavity_flow_fft_2d.py
