## idiot_proofer.py by Ryan Farber 20 December 2015
## Last modified: 20 December 2015 by Ryan Farber
"""
The purpose of this script is to check INPUTS* and other files
to make sure the inputs by the user are consistent (e.g., NY given
value of "NA" for 1D, etc.)

This currently assumes you're PWD is the problem directory.
To get the dimension (1,2,3) it assumes the problem directory name
follows the convention "ProblemName?D" where ? is 1,2,or 3.
"""
import os
import glob

DIM = os.getcwd().split("/")[-1][-2]
print DIM

file_to_read = glob.glob("INPUTS*")[0]

fid_read = open(file_to_read, "r")

#for line in fid_read:
  ## Do the error checking now
#  if 
