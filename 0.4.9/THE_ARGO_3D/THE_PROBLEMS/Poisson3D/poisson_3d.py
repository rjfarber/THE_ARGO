## poisson_3d.py
## by Ryan Farber 15 January 2015
## Last modified: 24 January 2015
"""
The purpose of this program is to apply the_argo to solve the
poisson equation for pressure in three dimensions.
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
ext = ".p"  # filename extension for pickling
p   = np.zeros( (NX,NY,NZ) ) # pressure
src = np.zeros( (NX,NY,NZ) ) # source term of poisson eqn

src[0.25*(NX-1), 0.25*(NY-1), : ] =  100
src[0.75*(NX-1), 0.75*(NY-1), : ] = -100


##Save the state of the initial condition of the fluid
cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, 'NA','NA','NA','NA',p,src,'NA'], file_name )


##Solve!
pct_err = 1.0; counter = 0; cycles += 1
while pct_err > my_fluid.MAE and counter < my_fluid.NI:
    ##prep for new iteration
    p_old = p.copy()
    counter += 1

    p = my_fluid.relax_pressure_poisson_3d(p,p_old, src)

    p = bs.set_bnds_fixed_3dX(p, 0, 0)
    p = bs.set_bnds_fixed_3dY(p, 0, 0)
    p = bs.set_bnds_fixed_3dZ(p, 0, 0)

    pct_err = (np.sum(np.abs(p) - np.abs(p_old))
            /  np.sum(np.abs(p_old)))
# end while
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)
save_state( [cycles, 'NA','NA','NA','NA',p,src,'NA'], file_name )


## end poisson_3d.py
