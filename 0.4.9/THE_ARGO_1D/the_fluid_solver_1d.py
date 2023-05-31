## the_fluid_solver_1d.py
## by Ryan Farber 20 January  2015
## Last modified: 20 December 2015

from the_fluid import The_Fluid

class The_Fluid_Solver_1D(The_Fluid):
    """
    The_Fluid_Solver_1D inherits __init__ from The_Fluid
    and requires minimally as input variables for instantiation:
    NX,NY,NZ (the number of zones in x,y, and z) and
    DX,DY,DZ (the width of a cell [zone] in x,y, and z).
    See an input file for further on instantiating a
    The_Fluid_Solver_1D instance.

    The_Fluid_Solver_1D includes methods for solving the navier stokes
    equations in three dimensions; however, the flows must be:
    constant density, constant viscosity, incompressible, and isothermal
    fluid flows.

    See the_bnd_setter_1d.py for setting boundary conditions and
    the_plotter_1d.py for plotting functions.
    """
    
    def linear_advect_1d(self, f):
        """
        Performs linear advection of a 1D field by backward differencing
        First Order.
        """

        f[1:] -= (self.DT/self.DX) * self.C * (f[1:] - f[:-1])
        
        return f
    # end linear_advect_1d

    def nonlinear_advect_1d(self, f):
        """
        Performs nonlinear advection of a 1D field by
        backward differencing. First Order
        """

        f[1:] -= (self.DT/self.DX) * f[1:] * (f[1:] - f[:-1])

        return f
    # end non_linear_advect_1d

    def diffuse_1d(self, f):
        """
        Performs diffusion of a 1D field by central differencing.
        Second Order
        """

        f[1:-1] += (self.DT/self.DX**2) * self.NU * (
                         f[2:] - 2*f[1:-1] + f[:-2])
        return f
    # end diffuse_1d
    
# end class The_Fluid_Solver_1D

## end the_fluid_solver_1d.py
