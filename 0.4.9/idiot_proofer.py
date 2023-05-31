## idiot_proofer.py by Ryan Farber 20 December 2015
## Last modified: 20 December 2015 by Ryan Farber
"""
The purpose of this script is to check INPUTS* and other files
to make sure the inputs by the user are consistent (e.g., NY given
value of "NA" for 1D, etc.)

This currently assumes the INPUT file has an import statement of the
form
    from the_fluid_solver_?d import The_Fluid_Solver_?D
where ? is 1,2, or 3
"""
import sys
import glob
import numpy as np

file_to_read = glob.glob("INPUTS*")[0]
fid_read = open(file_to_read, "r")

PLT_VARS="XMIN, XMAX, DATA_FILE_EXT, PLOT_FILE_EXT, LABEL, MAX_CYC_MAG"
PLT_VARS = ["my_fluid." + var in PLT_VARS]

print PLT_VARS
sys.exit()

DIM = 99
for line in fid_read:
  if "the_fluid_solver" in line:
    DIM = line[-2]
  # end if

  if "my_fluid.S" in line:
      if float(line.split()[2]) > 0.5:
        print("WARNING! You have set the stability constant to a value "
              "greater than 0.5; since the stability constant is the "
              "CFL number, the problem will blow up.")
      # end if
  elif DIM == "1":
    if np.intersect1d(["NY","NZ","DY","DZ"], line):
      if type(line.split()[2]) != str:
        print("WARNING! You have imported The_Fluid_Solver_1D but " +
              "NY,NZ,DY, or DZ is not of type string."
              "    It is advised you set each of those values to 'NA'"
              "OR import the correct dimension for your problem.")
      # end if
    elif np.intersect1d(PLT_VARS
    
    

  else:
    print("DEV NOTE: I haven't done idiot proofing for 2D or 3D yet.")
  # end if
# end for
    
## end idiot_proofer.py
