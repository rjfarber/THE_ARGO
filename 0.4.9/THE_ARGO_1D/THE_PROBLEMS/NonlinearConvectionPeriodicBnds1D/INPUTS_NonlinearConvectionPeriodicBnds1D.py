## INPUTS_NonlinearConvectionPeriodicBnds1D.py
## by Ryan Farber 21 January 2015
## Last modified: 22 January 2015
"""
This is an INPUTS file tailored for the NonlinearConvectionNoBnds1D
problem.
    NOTE: This file must be located two directories down from
the_fluid_solver_1d.py to function properly.
    ALSO NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_1d import The_Fluid_Solver_1D
from the_state_saver import save_state

##General Constatns
NX =  41   # number of zones in x
NY = "NA"  # not applicable to a 1D problem
NZ = "NA"

XMIN = 0.0
XMAX = 2.0

DX = (XMAX - XMIN) / (NX - 1)  # width of a cell in x
DY = "NA"
DZ = "NA"
##End General Constants


my_fluid = The_Fluid_Solver_1D(NX,NY,NZ, DX,DY,DZ)


##Problem Constants
my_fluid.NT   = 100        # number of time steps
my_fluid.DT   = DX/2.0     # chosen for stability
my_fluid.C    = 1.0        # the fluid wave speed
##End Problem Constants

my_fluid.SAVE_FREQ = 10    # save state every SAVE_FREQ cycles

##State Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX

my_fluid.LABEL       = "hat_periodic_boundaries"
my_fluid.MAX_CYC_MAG = len(str(my_fluid.NT))+1 # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
save_state(my_fluid, "my_fluid.p")


## end INPUTS_NonlinearConvectionPeriodicBnds1D.py
