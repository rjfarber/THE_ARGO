## INPUTS_Vortex3D.py
## by Ryan Farber 24 January 2015
## Last modified: 24 January 2015
"""
This is an INPUTS file tailored for the Vortex3D problem.

    NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_3d import The_Fluid_Solver_3D
from the_state_saver import save_state

##General Constatns
NX = 20   # number of zones in x
NY = 20   # number of zones in y
NZ = 20   # number of zones in z

XMIN = 0.0
XMAX = 1.0
YMIN = 0.0
YMAX = 1.0
ZMIN = 0.0
ZMAX = 1.0

DX = (XMAX - XMIN) / (NX)  # width of a cell in x
DY = (YMAX - YMIN) / (NY)  # width of a cell in y
DZ = (ZMAX - ZMIN) / (NZ)  # width of a cell in z
##End General Constants


my_fluid = The_Fluid_Solver_3D(NX,NY,NZ, DX,DY,DZ)


##Problem Constants
TF = 160  # time units simulation should complete
my_fluid.DT  = DX                    # chosen for stability
my_fluid.NT  = int(TF/my_fluid.DT)   # number of time steps
my_fluid.NI  =  200                  # number of iterations
my_fluid.RHO =  1.0                  # fluid density
##End Problem Constants


##Save Constants
my_fluid.SAVE_FREQ = 100       # save state every SAVE_FREQ cycles


##Plotting Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX
my_fluid.YMIN = YMIN
my_fluid.YMAX = YMAX
my_fluid.YMIN = ZMIN
my_fluid.YMAX = ZMAX

my_fluid.PLT_TYP     = "quiver and pcolormesh"
my_fluid.LABEL       = "vortex"
my_fluid.MAX_CYC_MAG = len(str(my_fluid.NT))+1 # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
save_state(my_fluid, "my_fluid.p")


## end INPUTS_Vortex3D.py
