## the_bnds_setter_1d.py
## by Ryan Farber 20 January 2015
## Last modified: 22 January 2015
"""
This file contains all the boundary condition functions to be used by
the_argo. Currently available boundary conditions are

- Fixed       (Dirichlet or Neumann or periodic, depending upon input)
- Periodic    (diffusion only)

and reflecting boundaries will hopefully one day be implemented.
"""

def set_bnds_fixed_1dI(f, amt):
    """Sets initial boundary cell of input 1d numpy array to a fixed
    amount of type int or float for Dirichlet (value) or Neumann (zero
    flux) or periodic boundary conditions."""
    
    f[ 0 ] = amt

    return f
# end set_bnds_fixed_1dI

def set_bnds_fixed_1dF(f, amt):
    """Sets final boundary cell of input 1d numpy array to a fixed
    amount of type int or float for Dirichlet (value) or Neumann (zero
    flux) or periodic boundary conditions."""
    
    f[ -1 ] = amt

    return f
# end set_bnds_fixed_1dF

def set_bnds_periodic_diffusion_1d(DT,DX,NU, f):
    """Sets periodic boundary conditions for diffusion."""

    f[ -1 ] += (DT/DX**2) * NU * ( f[ 0 ] - 2*f[ -1 ] + f[ -2 ] )
    
    return f
# end set_bnds_periodic_diffusion_1d

## the_bnds_setter_1d.py
