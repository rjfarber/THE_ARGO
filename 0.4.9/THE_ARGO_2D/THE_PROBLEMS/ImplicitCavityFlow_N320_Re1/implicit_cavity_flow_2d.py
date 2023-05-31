## implicit_cavity_flow_2d.py
## by Ryan Farber 15 February 2015
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
cycles = 0
steady_state = False
while not steady_state:
#for cycles in xrange(my_fluid.cycle_start, my_fluid.NT+1):
    u_old = u.copy()

    u[1:-1,1:-1] = my_fluid.nonlinear_advect_implicit_2d(u, u,v, XX,YY)
    v[1:-1,1:-1] = my_fluid.nonlinear_advect_implicit_2d(v, u,v, XX,YY)
    u[1:-1,1:-1] = my_fluid.diffuse_implicit_2d(u,u,my_fluid.NU)
    v[1:-1,1:-1] = my_fluid.diffuse_implicit_2d(v,v,my_fluid.NU)
    
    src[ 1:-1, 1:-1 ] = my_fluid.calc_source_2d(u,v)
    
    pct_err = 1.0; counter = 0
    while pct_err > my_fluid.MAE: #and counter < my_fluid.NI:
        ##prep for new iteration
        p_old = p.copy()
        counter += 1
        
        p = my_fluid.relax_pressure_poisson_2d(p, src)

        p[0,:] = p[1,:]; p[-1,:] = p[-2,:]; p[:,0] = p[:,1]; p[:,-1] = 0
        """
        p = my_fluid.set_bnds_fixed_2dXI( p, p[  1, : ] )
        p = my_fluid.set_bnds_fixed_2dXF( p, p[ -2, : ] )
        p = my_fluid.set_bnds_fixed_2dYI( p, p[  :, 1 ] )
        p = my_fluid.set_bnds_fixed_2dYF(p, 0.0)
        """

        ##calc percent error for this iteration
        pct_err = (np.sum(np.abs(p) - np.abs(p_old))
                /  np.sum(np.abs(p_old)))
    # end while

    u[1:-1,1:-1] -= my_fluid.apply_pressure_2dX(p, 1.0/my_fluid.RHO)
    v[1:-1,1:-1] -= my_fluid.apply_pressure_2dY(p, 1.0/my_fluid.RHO)
    
    ##Fix all boundary walls
    u[0,:] = 0; u[-1,:] = 0; u[:,0] = 0; u[:,-1] = 1
    v[0,:] = 0; v[-1,:] = 0; v[:,0] = 0; v[:,-1] = 0

    """    
    u = my_fluid.set_bnds_fixed_2d(u, 0)
    u = my_fluid.set_bnds_fixed_2dYF(u, 1)
    v = my_fluid.set_bnds_fixed_2d(v, 0)
    """

    cycles += 1

    pct_err_vel = (np.sum(np.abs(u) - np.abs(u_old))
                /  np.sum(np.abs(u_old)))

    if pct_err_vel < my_fluid.MAE:
        steady_state = True
    # end if

#    if cycles in [2, 200,700] or cycles % my_fluid.SAVE_FREQ == 0:
#        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
#                                           my_fluid.LABEL, ext)
#        save_state([cycles, u,v, "NA", p,src, "NA","NA", "NA","NA"],
#                   file_name ) 
    # end if
# end for
print("number of iterations = " + str(cycles) + " for Re = " +
      str(int(1.0/my_fluid.NU)) + " and N = " + str(my_fluid.NX))
file_name = my_fluid.get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                   my_fluid.LABEL, ext)
my_fluid.save_state([cycles, u,v, "NA", p,src, "NA","NA", "NA","NA"],
                   file_name )  


## end implicit_cavity_flow_2d.py
