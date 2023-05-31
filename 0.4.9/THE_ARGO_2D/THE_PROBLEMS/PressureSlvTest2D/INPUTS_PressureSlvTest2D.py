## INPUTS_PressureSlvTest2D.py
## by Ryan Farber 05 March 2015
## Last modified: 18 March 2015
"""
This is an INPUTS file tailored for the PressureSlvTest2D problem.

Uses only div_old to calculate pressure
Updates all zones in the mesh.

All values chosen to imitate:
https://perswww.kuleuven.be/~u0016541/MHD_sheets_pdf/AdvCFDles04.pdf

NOTE: NX,NY,NZ,DX,DY, and DZ are required constants to
instantiate the my_fluid object; additional problem specific
constants are added after creating the my_fluid object below.
"""
import numpy as np
import os; cur_path = os.getcwd()
import sys; sys.path.insert(0, '../..'); sys.path.insert(0, '../../..')
from the_fluid_solver_2d import The_Fluid_Solver_2D
from the_state_saver import save_state

##General Constants
NX =  51    # number of zones in x
NY =  51    # number of zones in y
NZ = "NA"    # not applicable to a 1D problem

XMIN = 0.0
XMAX = 2*np.pi
YMIN = 0.0
YMAX = 2*np.pi

DX = (XMAX - XMIN) / (NX-1)  # width of a cell in x
DY = (YMAX - YMIN) / (NY-1)  # width of a cell in y
DZ = "NA"
##End General Constants


my_fluid = The_Fluid_Solver_2D(NX,NY,NZ, DX,DY,DZ)


##Problem Constants
my_fluid.cycle_start = 1                # start cycle; use 1 as default
my_fluid.S   = 1E-1                     # stability constant (vmax < 5E2)
my_fluid.DT  = DX*my_fluid.S            # chosen for stability
my_fluid.NT  = 1600                     # number of time steps
my_fluid.NI  = int(2e2)                 # number of iterations

my_fluid.RHO = 1.0                      # density of the fluid
my_fluid.MU  = 1.0                      # magnetic permeability
my_fluid.NU  = 0.04                     # viscous diffusion
my_fluid.ETA = 0.04                     # resistivity (mag diffusion)
##End Problem Constants


my_fluid.SAVE_SP   = range(100) + range(100,201,10) # save first x
my_fluid.SAVE_FREQ = 100      # save state every SAVE_FREQ cycles


##Plotting Constants
my_fluid.XMIN = XMIN
my_fluid.XMAX = XMAX
my_fluid.YMIN = YMIN
my_fluid.YMAX = YMAX

my_fluid.PLT_TYP     = "quiver pcolormesh concentration and magnetism"
my_fluid.LABEL       = "hat"
my_fluid.MAX_CYC_MAG = 9                     # magnitude of max cycles
##End Plotting Constants


os.chdir(cur_path + "/StateFiles")
save_state(my_fluid, "my_fluid.p")


## end INPUTS_PressureSlvTest2D.py
