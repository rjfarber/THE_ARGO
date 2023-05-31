## INPUTS_JacobiVsFFT2D.py
## by Ryan Farber 23 January 2015
## Last modified: 01 April   2015
"""
This is an INPUTS file tailored for the JacobiVsFFT2D problem.

NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_2d import The_Fluid_Solver_2D
from the_state_saver import save_state

##General Constatns
NX =  20   # number of zones in x
NY =  20   # number of zones in y
NZ = "NA"  # not applicable to a 1D problem

XMIN = 0.0
XMAX = 2.0
YMIN = 0.0
YMAX = 2.0

DX = (XMAX - XMIN) / (NX - 2)  # width of a cell in x
DY = (YMAX - YMIN) / (NY - 2)  # width of a cell in y
DZ = "NA"
##End General Constants


my_fluid = The_Fluid_Solver_2D(NX,NY,NZ, DX,DY,DZ)


##Problem Constants
my_fluid.cycle_start = 1      # start cycle; use 1 for default
my_fluid.NI  = 100            # number of iterations
my_fluid.MAE = 1E-2           # maximum allowed error
##End Problem Constants


##Plotting Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX
my_fluid.YMIN = YMIN
my_fluid.YMAX = YMAX

my_fluid.PLT_TYP     = "surface pressure exclude boundaries"
my_fluid.LABEL       = "pressure_test"
my_fluid.MAX_CYC_MAG = 2                     # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
save_state(my_fluid, "my_fluid.p")


## end INPUTS_JacobiVsFFT2D.py
