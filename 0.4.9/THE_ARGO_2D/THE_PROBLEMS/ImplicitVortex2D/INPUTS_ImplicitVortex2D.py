## INPUTS_ImplicitVortex2D.py
## by Ryan Farber 23 January 2015
## Last modified: 02 April   2015
"""
This is an INPUTS file tailored for the ImplicitVortex2D problem.

NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_2d import The_Fluid_Solver_2D
from the_state_saver import save_state

##General Constatns
NX =  21   # number of zones in x
NY =  21   # number of zones in y
NZ = "NA"  # not applicable to a 1D problem

XMIN = 0.0
XMAX = 1.0
YMIN = 0.0
YMAX = 1.0

DX = (XMAX - XMIN) / (NX - 1)  # width of a cell in x
DY = (YMAX - YMIN) / (NY - 1)  # width of a cell in y
DZ = "NA"
##End General Constants


my_fluid = The_Fluid_Solver_2D(NX,NY,NZ, DX,DY,DZ)


##Problem Constants
my_fluid.cycle_start = 1                # start cycle; use 1 as default
my_fluid.S   = 0.2                      # stability constant
my_fluid.DT  = DX*my_fluid.S            # chosen for stability
my_fluid.NT  = 3200                     # number of time steps
my_fluid.NI  = 500                      # number of iterations

my_fluid.NU  = 0.2                      # viscous diffusion
my_fluid.RHO = 1.0                      # density of the fluid
##End Problem Constants


my_fluid.SAVE_SP   = range(10) + range(10,101,10) # save first x
my_fluid.SAVE_FREQ = 100      # save state every SAVE_FREQ cycles


##Plotting Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX
my_fluid.YMIN = YMIN
my_fluid.YMAX = YMAX

my_fluid.PLT_TYP     = "quiver and pcolormesh concentration"
my_fluid.LABEL       = "vortex"
my_fluid.MAX_CYC_MAG = 9                     # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
save_state(my_fluid, "my_fluid.p")


## end INPUTS_ImplicitVortex2D.py
