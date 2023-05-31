## explicit_cavity_flow_2d.py
## by Ryan Farber 30 December 2014
## Last modified: 05 May      2015
"""
The purpose of this program is to apply the_argo to solve the
navier stokes equations in two dimensions for cavity flow.
"""
import glob
import numpy as np
import sys; sys.path.insert(0, "../../..")    
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Setup
ext = ".p"  # filename extension for pickling
BETA = 1.8

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
    u = the_data[1]; v = the_data[2]
    the_data = 0 # to save memory
# end if


##Solve!
for cycles in xrange(my_fluid.cycle_start, my_fluid.NT+1):
      ## First advect and diffuse velocity field
    u_old = u.copy(); v_old = v.copy()
    
    u[1:-1,1:-1] = u_old[1:-1,1:-1] - my_fluid.advect_explicit_CFD_2dX(u_old,v_old) \
                 + my_fluid.diffuse_CFD_2d(u_old)
    v[1:-1,1:-1] = v_old[1:-1,1:-1] + my_fluid.advect_explicit_CFD_2dY(u_old,v_old) \
                 + my_fluid.diffuse_CFD_2d(v_old)

    ##Fix all boundary walls
    u[0,:] = 0; u[-1,:] = 0; u[:,0] = 0; u[:,-1] = 1
    v[0,:] = 0; v[-1,:] = 0; v[:,0] = 0; v[:,-1] = 0

        ## Solid wall BCs
      ## XI
#    u[ 0,:] = -u[ 1, :]      v[ 0,:] =  v[ 1,:]
      ## XF
#    u[-1,:] = -u[-2, :]      v[-1,:] =  v[-2,:]
      ## YI
#    u[ :,0] =  u[ :, 1]      v[ :,0] = -v[ :,1]
      ## YF
#    u[:,-1] =  u[ :,-2]      v[:,-1] = -v[:,-2]
    
      ## Second, compute the pressure
    pct_err = 1.0; counter = 0
    while counter < 25:
        ##prep for new iteration
        p_old = p.copy()
        counter += 1
        
        p = my_fluid.SOR_pressure_CFD_2d(p, u,v, BETA)

        ##calc percent error for this iteration
        pct_err = (np.sum(np.abs(p) - np.abs(p_old))
                /  np.sum(np.abs(p_old)))
    # end while
    print("pct_err: ", pct_err)
          
        ## Third, pressure-correct the velocity (make div(vec{u}) = 0)

    u[1:-1,1:-1] -= (my_fluid.DT/my_fluid.DX)*(p[2:,1:-1] - p[1:-1,1:-1])
    v[1:-1,1:-1] -= (my_fluid.DT/my_fluid.DY)*(p[1:-1,2:] - p[1:-1,1:-1])

    div = 0.5*(u[2:,1:-1] - u[:-2,1:-1])/my_fluid.DX + \
          0.5*(v[1:-1,2:] - v[1:-1,:-2])/my_fluid.DY

    print("divergence of velocity: ", np.sum(np.abs(div)))

    raw_input()

        ## Solid wall BCs
      ## XI
#    u[ 0,:] = -u[ 1, :];      v[ 0,:] =  v[ 1,:]
      ## XF
#    u[-1,:] = -u[-2, :];      v[-1,:] =  v[-2,:]
      ## YI
#    u[ :,0] =  u[ :, 1];      v[ :,0] = -v[ :,1]
      ## YF
#    u[:,-1] =  1.0;           v[:,-1] = -v[:,-2]
    
        
    ##Fix all boundary walls
    u[0,:] = 0; u[-1,:] = 0; u[:,0] = 0; u[:,-1] = 1
    v[0,:] = 0; v[-1,:] = 0; v[:,0] = 0; v[:,-1] = 0
    
#    u = my_fluid.set_bnds_fixed_2d(u, 0)
#    u = my_fluid.set_bnds_fixed_2dYF(u, 1)
#    v = my_fluid.set_bnds_fixed_2d(v, 0)   

    if cycles in [2, 200,700] or cycles % my_fluid.SAVE_FREQ == 0:
        file_name = my_fluid.get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                           my_fluid.LABEL, ext)
        my_fluid.save_state([cycles, u,v, "NA", p,src, "NA","NA", "NA","NA"],
                   file_name ) 
    # end if
# end for
pct_err = (np.sum(np.abs(u) - np.abs(u_old))
                /  np.sum(np.abs(u_old)))
print("pct_err vel: ", pct_err)

file_name = my_fluid.get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                   my_fluid.LABEL, ext)
my_fluid.save_state([cycles, u,v, "NA", p,src, "NA","NA", "NA","NA"],
           file_name ) 

## end explicit_cavity_flow_2d.py
