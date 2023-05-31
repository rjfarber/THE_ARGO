## INPUTS_ExplicitBurgers2D.py
## by Ryan Farber 23 January 2015
## Last modified: 07 March   2015
"""
This is an INPUTS file tailored for the ExplicitBurgers2D problem.

NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_2d import The_Fluid_Solver_2D

##General Constatns
NX =  41   # number of zones in x
NY =  41   # number of zones in y
NZ = "NA"  # not applicable to a 1D problem

XMIN = 0.0
XMAX = 2.0
YMIN = 0.0
YMAX = 2.0

DX = (XMAX - XMIN) / (NX - 1)  # width of a cell in x
DY = (YMAX - YMIN) / (NY - 1)  # width of a cell in y
DZ = "NA"
##End General Constants


my_fluid = The_Fluid_Solver_2D(NX,NY,NZ, DX,DY,DZ)


##Problem Constants
my_fluid.cycle_start = 1                        # start cycle; use 1
my_fluid.NT = 9600                              # number of time steps
my_fluid.NU = 0.01                              # viscous diffusion
my_fluid.S  = 9E-4                              # stability constant
my_fluid.DT = (my_fluid.S*DX*DY) / my_fluid.NU  # chosen for stability
##End Problem Constants


my_fluid.SAVE_FREQ = 600      # save state every SAVE_FREQ cycles


##Plotting Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX
my_fluid.YMIN = YMIN
my_fluid.YMAX = YMAX

my_fluid.PLT_TYP     = "wireframe velocity"
my_fluid.LABEL       = "diffuse_hat_fixed_boundaries"
my_fluid.MAX_CYC_MAG = 9                     # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
my_fluid.save_state(my_fluid, "my_fluid.p")


## end INPUTS_ExplicitBurgers2D.py
