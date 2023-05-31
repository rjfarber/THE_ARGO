## INPUTS_Poisson3D.py
## by Ryan Farber 24 January 2015
## Last modified: 24 January 2015
"""
This is an INPUTS file tailored for the Poisson3D problem.

    NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_3d import The_Fluid_Solver_3D
from the_state_saver import save_state

##General Constatns
NX = 51   # number of zones in x
NY = 51   # number of zones in y
NZ = 51   # number of zones in z

XMIN = 0.0
XMAX = 2.0
YMIN = 0.0
YMAX = 2.0
ZMIN = 0.0
ZMAX = 2.0

DX = (XMAX - XMIN) / (NX - 1)  # width of a cell in x
DY = (YMAX - YMIN) / (NY - 1)  # width of a cell in y
DZ = (ZMAX - ZMIN) / (NZ - 1)  # width of a cell in z
##End General Constants


my_fluid = The_Fluid_Solver_3D(NX,NY,NZ, DX,DY,DZ)


##Problem Constants
my_fluid.NI  =  100   # number of iterations
my_fluid.MAE = 1E-2   # maximum allowed error
##End Problem Constants


##Plotting Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX
my_fluid.YMIN = YMIN
my_fluid.YMAX = YMAX
my_fluid.YMIN = ZMIN
my_fluid.YMAX = ZMAX

my_fluid.PLT_TYP     = "surface pressure"
my_fluid.LABEL       = "poisson"
my_fluid.MAX_CYC_MAG = 2                  # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
save_state(my_fluid, "my_fluid.p")


## end INPUTS_Poisson3D.py