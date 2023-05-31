## laplace_3d.py
## by Ryan Farber 15 January 2015
## Last modified: 24 January 2015
"""
The purpose of this program is to apply the_argo to solve the
laplace equation for pressure in three dimensions.
"""
import numpy as np
import sys; sys.path.insert(0, "../../");import the_bnds_setter_3d as bs 
sys.path.insert(0, "../../..")    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))
NX = my_fluid.NX; NY = my_fluid.NY; NZ = my_fluid.NZ

##Setup
Y = np.linspace(my_fluid.YMIN, my_fluid.YMAX, NY)
ext = ".p"  # filename extension for pickling
p   = np.zeros( (NX,NY,NZ) ) # pressure
src = np.zeros( (NX,NY,NZ) ) # source term of poisson eqn

p[-1,:,:] = Y


##Save the state of the initial condition of the fluid
cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, 'NA','NA','NA','NA',p,src,'NA'], file_name )


##Solve!
pct_err = 1.0
while pct_err > my_fluid.MAE:
    p_old = p.copy()

    ##Solves laplace eqn for src=0
    p = my_fluid.relax_pressure_poisson_3d(p,p_old, src)

    ##Boundary conditions
        ##p = 0 @ x = 0 and p = y at x = 2
    p[0,:,:] = 0; p[-1,:,:] = Y
        ##dp/dy = 0 at y = 0,2 (partial derivs)
    p[:,0,:] = p[:,1,:]; p[:,-1,:] = p[:,-2,:]

    pct_err = (np.sum(np.abs(p) - np.abs(p_old))
            /  np.sum(np.abs(p_old)))
# end while
cycles = 2; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, 'NA','NA','NA','NA',p,src,'NA'], file_name )


## end laplace_3d.py
