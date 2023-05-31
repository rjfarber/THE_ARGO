## INPUTS_DoubleVortexOppositeSign2D.py
## by Ryan Farber 27 January 2015
## Last modified: 27 January 2015
"""
This is an INPUTS file tailored for the DoubleVortexOppositeSign2D
problem.

VortexA uses only div_old to calculate pressure
VortexA updates all zones in the mesh.

    NOTE: This file must be located two directories down from
the_fluid_solver_2d.py to function properly.

    ALSO NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_2d import The_Fluid_Solver_2D
from the_state_saver import save_state

##General Constatns
NX =  51   # number of zones in x
NY =  51   # number of zones in y
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
TF = 80*160                             # time units to complete
my_fluid.S   = 10.0                     # stability constant
my_fluid.DT  = my_fluid.S*DX            # chosen for stability
my_fluid.NT  = int(TF/my_fluid.DT)      # number of time steps
my_fluid.NI  = 200                      # number of iterations
my_fluid.RHO = 1.0                      # density of the fluid
##End Problem Constants


##Save Constants
my_fluid.SAVE_FREQ = 500      # save state every SAVE_FREQ cycles


##Plotting Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX
my_fluid.YMIN = YMIN
my_fluid.YMAX = YMAX

my_fluid.PLT_TYP     = "quiver and pcolormesh"
my_fluid.LABEL       = "double_vortex_opposite_sign_S=10.0"
my_fluid.MAX_CYC_MAG = len(str(my_fluid.NT))+1 # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
save_state(my_fluid, "my_fluid.p")


## end INPUTS_DoubleVortexOppositeSign2D.py
