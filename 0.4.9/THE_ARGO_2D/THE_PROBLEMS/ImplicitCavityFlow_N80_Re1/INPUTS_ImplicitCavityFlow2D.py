## INPUTS_ImplicitCavityFlow2D.py
## by Ryan Farber 15 February 2015
## Last modified: 05 May      2015
"""
This is an INPUTS file tailored for the ImplicitCavityFlow2D problem.

NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_2d import The_Fluid_Solver_2D

##General Constatns
NX =  81   # number of zones in x (plus one)
NY =  81   # number of zones in y (plus one)
NZ = "NA"   # not applicable to a 1D problem

XMIN = 0.0
XMAX = 2.0
YMIN = 0.0
YMAX = 2.0

DX = (XMAX - XMIN) / (NX - 1)  # width of a cell in x
DY = (YMAX - YMIN) / (NY - 1)  # width of a cell in y
DZ = "NA"
##End General Constants


my_fluid = The_Fluid_Solver_2D(NX,NY,NZ, DX,DY,DZ)


Re = 1

##Problem Constants
my_fluid.cycle_start = 1      # start cycle; use 1 for default
my_fluid.NT  = 800            # number of time steps
my_fluid.NI  = 50             # number of iterations
my_fluid.DT  = 1E-3           
my_fluid.NU  = 1.0/Re         # viscous diffusion
my_fluid.RHO = 1.0            # density of the fluid
my_fluid.MAE = 1E-4           # maximum allowed error
##End Problem Constants

my_fluid.SAVE_FREQ = 8

##Plotting Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX
my_fluid.YMIN = YMIN
my_fluid.YMAX = YMAX

my_fluid.PLT_TYP     = "quiver contours include boundaries"
my_fluid.LABEL       = "implicit_cavity_flow_N" + str(NX) + "_Re" + str(Re)
my_fluid.MAX_CYC_MAG = 9                     # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
my_fluid.save_state(my_fluid, "my_fluid.p")


## end INPUTS_ImplicitCavityFlow2D.py
