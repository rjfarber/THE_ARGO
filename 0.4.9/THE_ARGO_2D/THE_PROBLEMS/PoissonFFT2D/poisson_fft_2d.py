## poisson_2d.py
## by Ryan Farber 30 December 2014
## Last modified: 18 March    2015
"""
The purpose of this program is to apply the_argo to solve the
poisson equation for pressure in two dimensions.
"""
import numpy as np
import sys; sys.path.insert(0, "../../..")    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Setup
ext = ".p"  # filename extension for pickling
p   = np.zeros( (my_fluid.NX,my_fluid.NY) ) # pressure
src = np.zeros( (my_fluid.NX,my_fluid.NY) ) # source term of poisson eqn

src[0.25*(my_fluid.NX-1), 0.25*(my_fluid.NY-1)] =  100
src[0.75*(my_fluid.NX-1), 0.75*(my_fluid.NY-1)] = -100


##Save the state of the initial condition of the fluid
cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                               my_fluid.LABEL, ext)
save_state([cycles, "NA", "NA","NA", p,src, "NA","NA", "NA","NA"],
           file_name )

##Solve!
pct_err = 1.0; counter = 0; cycles += 1
while pct_err > my_fluid.MAE and counter < my_fluid.NI:
    ##prep for new iteration
    p_old = p.copy()
    counter += 1

    p = my_fluid.transform_pressure_poisson_2d(p, src)
    
    ##Fix all boundary walls at 0.0
    p = my_fluid.set_bnds_fixed_2d(p, 0)

    pct_err = (np.sum(np.abs(p) - np.abs(p_old))
            /  np.sum(np.abs(p_old)))
# end while
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                   my_fluid.LABEL, ext)
save_state([cycles, "NA", "NA","NA", p,src, "NA","NA", "NA","NA"],
           file_name )

## end poisson_2d.py
