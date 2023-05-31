## INPUTS_Burgers1D.py
## by Ryan Farber 21 January 2015
## Last modified: 22 January 2015
"""
This is an INPUTS file tailored for the Diffusion1D problem.
    NOTE: This file must be located two directories down from
the_fluid_solver_1d.py to function properly.
    ALSO NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import numpy as np
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_1d import The_Fluid_Solver_1D
from the_state_saver import save_state

##General Constatns
NX =  101  # number of zones in x
NY = "NA"  # not applicable to a 1D problem
NZ = "NA"

XMIN = 0.0
XMAX = 2*np.pi

DX = (XMAX - XMIN) / (NX - 1)  # width of a cell in x
DY = "NA"
DZ = "NA"
##End General Constants


my_fluid = The_Fluid_Solver_1D(NX,NY,NZ, DX,DY,DZ)


##Problem Constants
my_fluid.NT = 100                            # number of time steps
my_fluid.NU = 0.07                           # viscous diffusion
my_fluid.DT = my_fluid.DX*my_fluid.NU        # chosen for stability
##End Problem Constants


##State Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX

my_fluid.LABEL       = "sawtooth"
my_fluid.MAX_CYC_MAG = len(str(my_fluid.NT))+1 # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
save_state(my_fluid, "my_fluid.p")


## end INPUTS_Burgers1D.py
