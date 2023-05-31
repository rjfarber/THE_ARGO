## the_fluid_solver_3d.py
## by Ryan Farber 20 January 2015
## Last modified: 27 January 2015

from the_fluid import The_Fluid

class The_Fluid_Solver_3D(The_Fluid):
    """
    The_Fluid_Solver_3D inherits __init__ from The_Fluid
    and requires minimally as input variables for instantiation:
    NX,NY,NZ (the number of zones in x,y, and z) and
    DX,DY,DZ (the width of a cell [zone] in x,y, and z).
    See an input file for further on instantiating a
    The_fluid_Solver_3D instance.

    The_Fluid_Solver_3D includes methods for solving the navier stokes
    equations in three dimensions; however, the flows must be:
    constant density, constant viscosity, incompressible, and isothermal
    fluid flows.

    See the_bnds_setter.py for setting boundaries, see the_plotter_3d.py
    for plotting and the_div_calcer.py for calculating the divergence.
    """
    
##3D SOLVER FUNCTIONS
    def linear_convect_3d(self, f,f_old):
        """Performs linear convection of a 3D field by
        backward differencing."""
            ##first in x
        f[ 1:, 1:, 1: ] -= (self.DT/self.DX)*self.C*(
                f_old[ 1:, 1:, 1: ] - f_old[  :-1, 1:  , 1:   ])

            ##now in y        
        f[ 1:, 1:, 1: ] -= (self.DT/self.DY)*self.C*(
                f_old[ 1:, 1:, 1: ] - f_old[ 1:  ,  :-1, 1:   ])

            ##last in z       
        f[ 1: ,1: ,1: ] -= (self.DT/self.DZ)*self.C*(
                f_old[ 1:, 1:, 1: ] - f_old[ 1:  , 1:  ,  :-1 ])
        return f
    # end linear_convect_3d

    def nonlinear_convect_3dX(self, fx, fx_old,fy_old,fz_old):
        """Performs nonlinear convection of the x-component of
        a 3D field by backward differencing."""
            ##first in x
        fx[ 1:, 1:, 1: ] -= (fx_old[ 1:, 1:, 1: ]
 * (self.DT/self.DX) * (fx_old[ 1:  , 1:  , 1:   ]
                     -  fx_old[  :-1, 1:  , 1:   ]))
        
            ##now in y
        fx[ 1:, 1:, 1: ] -= (fy_old[ 1:, 1:, 1: ]
 * (self.DT/self.DY) * (fx_old[ 1:  , 1:  , 1:   ]
                     -  fx_old[ 1:  ,  :-1, 1:   ]))

            ##last in z
        fx[ 1: ,1: ,1: ] -= (fz_old[ 1:, 1:, 1: ]
 * (self.DT/self.DZ) * (fx_old[ 1:  , 1:  , 1:   ]
                     -  fx_old[ 1:  , 1:  ,  :-1 ]))
        return fx
    # end nonlinear_convect_3dX

    def nonlinear_convect_3dY(self, fy, fx_old,fy_old,fz_old):
        """Performs nonlinear convection of the y-component of
        a 3D field by backward differencing."""
            ##first in x
        fy[ 1:, 1:, 1: ] -= (fx_old[ 1:, 1:, 1: ]
 * (self.DT/self.DX) * (fy_old[ 1:  , 1:  , 1:   ]
                     -  fy_old[  :-1, 1:  , 1:   ]))
        
            ##now in y
        fy[ 1:, 1:, 1: ] -= (fy_old[ 1:, 1:, 1: ]
 * (self.DT/self.DY) * (fy_old[ 1:  , 1:  , 1:   ]
                     -  fy_old[ 1:  ,  :-1, 1:   ]))

            ##last in z
        fy[ 1: ,1: ,1: ] -= (fz_old[ 1:, 1:, 1: ]
 * (self.DT/self.DZ) * (fy_old[ 1:  , 1:  , 1:   ]
                     -  fy_old[ 1:  , 1:  ,  :-1 ]))
        return fy
    # end nonlinear_convect_3dY

    def nonlinear_convect_3dZ(self, fz, fx_old,fy_old,fz_old):
        """Performs nonlinear convection of the z-component of
        a 3D field by backward differencing."""
            ##first in x
        fz[ 1:, 1:, 1: ] -= (fx_old[ 1:, 1:, 1: ]
 * (self.DT/self.DX) * (fz_old[ 1:  , 1:  , 1:   ]
                     -  fz_old[  :-1, 1:  , 1:   ]))
        
            ##now in y
        fz[ 1:, 1:, 1: ] -= (fy_old[ 1:, 1:, 1: ]
 * (self.DT/self.DY) * (fz_old[ 1:  , 1:  , 1:   ]
                     -  fz_old[ 1:  ,  :-1, 1:   ]))

            ##last in z
        fz[ 1: ,1: ,1: ] -= (fz_old[ 1:, 1:, 1: ]
 * (self.DT/self.DZ) * (fz_old[ 1:  , 1:  , 1:   ]
                     -  fz_old[ 1:  , 1:  ,  :-1 ]))
        return fz
    # end nonlinear_convect_3dZ
    
    def diffuse_3d(self, f,f_old):
        """Performs diffusion of a 3D field by central differencing;
        viscosity is assumed to be constant."""
            ##first in x
        f[ 1:-1, 1:-1, 1:-1 ] += self.NU*(self.DT/self.DX**2
 * (f_old[ 2:  , 1:-1, 1:-1 ] - 2*f_old[ 1:-1, 1:-1, 1:-1 ]
 +  f_old[  :-2, 1:-1, 1:-1 ]))

            ##now in y
        f[ 1:-1, 1:-1, 1:-1 ] += self.NU*(self.DT/self.DY**2
 * (f_old[ 1:-1, 2:  , 1:-1 ] - 2*f_old[ 1:-1, 1:-1, 1:-1 ]
 +  f_old[ 1:-1,  :-2, 1:-1 ]))

            ##now in z
        f[ 1:-1, 1:-1, 1:-1 ] += self.NU*(self.DT/self.DZ**2
 * (f_old[ 1:-1, 1:-1, 2:   ] - 2*f_old[ 1:-1, 1:-1, 1:-1 ]
 +  f_old[ 1:-1, 1:-1,  :-2 ]))

        return f
    # end diffuse_3d

    def relax_diffusion_3d(self, f):
        """Performs diffusion of a 3D field by implicit
        central differencing; viscosity is assumed to be constant.
            NOTE: I can maybe just use relax_pressure_poisson if
        I can generalize it.
            WHY: is it A_ETA / (1+6*A_ETA)? Shouldn't it be just
        divided by 6? (see Implicit Diffusion page of notebook)"""

        raise("error: not operational")

        f[ 1:-1, 1:-1, 1:-1 ] += A_ETA / (1+6.0*A_ETA)*(

                                f[ 2:  , 1:-1, 1:-1 ]
                              + f[  :-2, 1:-1, 1:-1 ]
                              + f[ 1:-1, 2:  , 1:-1 ]
                              + f[ 1:-1,  :-2, 1:-1 ]
                              + f[ 1:-1, 1:-1, 2:   ]
                              + f[ 1:-1, 1:-1,  :-2 ]
                                                       )
        return f

    def advect_3d(self, f, i0,i1, j0,j1, k0,k1, s0,s1, t0,t1, u0,u1):
        """Advects a 3D (velocity) field by cell-centered back-tracking
        and applying necessary weights."""
        
        f[1:-1,1:-1,1:-1] = (
                              s0 * (
                                     t0 * 
               ( u0 * f[ i0, j0, k0 ] + u1 * f[ i0, j0, k1 ] )
                                   +  
                                     t1 *
               ( u0 * f[ i0, j1, k0 ] + u1 * f[ i0, j1, k1 ] )
                                   )
                                   +
                              s1 * (
                                     t0 *
               ( u0 * f[ i1, j0 ,k0 ] + u1 * f[ i1, j0, k1 ] )
                                   +
                                     t1 *
               ( u0 * f[ i1, j1, k0 ] + u1 * f[ i1, j1, k1 ] )
                                   )
                            )
        return f
    # end advect_3d

    def apply_pressure_3dX(self, fx, p):
        """applies the x-component of the pressure gradient to the
        x-component of a 3D field [velocity] by central differencing;
        density is assumed to be constant."""

        fx[ 1:-1, 1:-1, 1:-1 ] -= self.DT/(2*self.RHO*self.DX)*(
            p[ 2:  , 1:-1, 1:-1 ] - p[  :-2, 1:-1, 1:-1 ])
        
        return fx
    # end apply_pressure_3dX

    def apply_pressure_3dY(self, fy, p):
        """applies the y-component of the pressure gradient to the
        y-component of a 3D field [velocity] by central differencing;
        density is assumed to be constant."""
        
        fy[ 1:-1, 1:-1, 1:-1 ] -= self.DT/(2*self.RHO*self.DY)*(
            p[ 1:-1, 2:  , 1:-1 ] - p[ 1:-1,  :-2, 1:-1 ])

        return fy
    # end apply_pressure_3dY

    def apply_pressure_3dZ(self, fz, p):
        """applies the z-component of the pressure gradient to the
        z-component of a 3D field [velocity] by central differencing;
        density is assumed to be constant."""
        
        fz[ 1:-1, 1:-1, 1:-1 ] -= self.DT/(2*self.RHO*self.DZ)*(
            p[ 1:-1, 1:-1, 2:   ] - p[ 1:-1, 1:-1,  :-2 ])

        return fz
    # end apply_pressure_3dZ

    def apply_force_3d(self, f, g):
        """applies a component of the acceleration due to a force [such
        as gravity] (g) to a component of a 3D field [velocity] (f)."""
        
        f[ 1:-1, 1:-1, 1:-1 ] += self.DT * g[ 1:-1, 1:-1, 1:-1 ]

        return f
    # end apply_force_3d

    def calc_source_3d(self, src, u,v,w):
        """Calculates source for 3D poisson equation
        of navier-stokes."""

        src[1:-1, 1:-1, 1:-1] = self.RHO*(
1.0/self.DT*(
                    (u[ 2:  , 1:-1, 1:-1 ] - u[  :-2, 1:-1, 1:-1 ])
                                            / (2*self.DX)  +
                    (v[ 1:-1, 2:  , 1:-1 ] - v[ 1:-1,  :-2, 1:-1 ])
                                            / (2*self.DY)  +
                    (w[ 1:-1, 1:-1, 2:   ] - w[ 1:-1, 1:-1,  :-2 ])
                                            / (2*self.DZ)
            )
              -
                (   (u[ 2:  , 1:-1, 1:-1 ] - u[  :-2, 1:-1, 1:-1 ])
                                            / (2*self.DX)
                )**2
              -
                (   (v[ 1:-1, 2:  , 1:-1 ] - v[ 1:-1,  :-2, 1:-1 ])
                                            / (2*self.DY)
                )**2
              -
                (   (w[ 1:-1, 1:-1, 2:   ] - w[ 1:-1, 1:-1,  :-2 ])
                                            / (2*self.DZ)
                )**2
              -
              2*(   (v[ 2:  , 1:-1, 1:-1 ] - v[  :-2, 1:-1, 1:-1 ])
                                            / (2*self.DX)
                  * (u[ 1:-1, 2:  , 1:-1 ] - u[ 1:-1,  :-2, 1:-1 ])
                                            / (2*self.DY)
                )
              -
              2*(   (w[ 2:  , 1:-1, 1:-1 ] - w[  :-2, 1:-1, 1:-1 ])
                                            / (2*self.DX)
                  * (u[ 1:-1, 1:-1, 2:   ] - u[ 1:-1, 1:-1,  :-2 ])
                                            / (2*self.DZ)
                )
              -
              2*(   (w[ 1:-1, 2:  , 1:-1 ] - w[ 1:-1,  :-2, 1:-1 ])
                                            / (2*self.DY)
                  * (v[ 1:-1, 1:-1, 2:   ] - v[ 1:-1, 1:-1,  :-2 ])
                                            / (2*self.DZ)
                )
                                          )
        return src
    # end calc_source_3d

    def calc_source_div_only_3d(self, src, u,v,w):
        """Calculates source for 3D poisson equation
        of navier-stokes, including only the divergence term."""

        src[1:-1, 1:-1, 1:-1] = self.RHO*(
1.0/self.DT*(
                    (u[ 2:  , 1:-1, 1:-1 ] - u[  :-2, 1:-1, 1:-1 ])
                                            / (2*self.DX)  +
                    (v[ 1:-1, 2:  , 1:-1 ] - v[ 1:-1,  :-2, 1:-1 ])
                                            / (2*self.DY)  +
                    (w[ 1:-1, 1:-1, 2:   ] - w[ 1:-1, 1:-1,  :-2 ])
                                            / (2*self.DZ)
            )
                                          )
        return src
    # end calc_source_div_only_3d

    def calc_source_all_terms_3d(self, src, u,v,w):
        """Calculates source for all terms of the 3D poisson equation
        of navier-stokes; NOT OPERATIONAL"""

        raise("error: not operational")

        src[1:-1, 1:-1, 1:-1] = self.RHO*(
1.0/self.DT*(
                    (u[ 2:  , 1:-1, 1:-1 ] - u[  :-2, 1:-1, 1:-1 ])
                                            / (2*self.DX)  +
                    (v[ 1:-1, 2:  , 1:-1 ] - v[ 1:-1,  :-2, 1:-1 ])
                                            / (2*self.DY)  +
                    (w[ 1:-1, 1:-1, 2:   ] - w[ 1:-1, 1:-1,  :-2 ])
                                            / (2*self.DZ)
            )
              -
                (   (u[ 2:  , 1:-1, 1:-1 ] - u[  :-2, 1:-1, 1:-1 ])
                                            / (2*self.DX)
                )**2
              -
                (   (v[ 1:-1, 2:  , 1:-1 ] - v[ 1:-1,  :-2, 1:-1 ])
                                            / (2*self.DY)
                )**2
              -
                (   (w[ 1:-1, 1:-1, 2:   ] - w[ 1:-1, 1:-1,  :-2 ])
                                            / (2*self.DZ)
                )**2
              -
              2*(   (v[ 2:  , 1:-1, 1:-1 ] - v[  :-2, 1:-1, 1:-1 ])
                                            / (2*self.DX)
                  * (u[ 1:-1, 2:  , 1:-1 ] - u[ 1:-1,  :-2, 1:-1 ])
                                            / (2*self.DY)
                )
              -
              2*(   (w[ 2:  , 1:-1, 1:-1 ] - w[  :-2, 1:-1, 1:-1 ])
                                            / (2*self.DX)
                  * (u[ 1:-1, 1:-1, 2:   ] - u[ 1:-1, 1:-1,  :-2 ])
                                            / (2*self.DZ)
                )
              -
              2*(   (w[ 1:-1, 2:  , 1:-1 ] - w[ 1:-1,  :-2, 1:-1 ])
                                            / (2*self.DY)
                  * (v[ 1:-1, 1:-1, 2:   ] - v[ 1:-1, 1:-1,  :-2 ])
                                            / (2*self.DZ)
                )
##The below assumes DX = DY = DZ; also I'm not sure how to adjust
##the above to be consistent with the below.
+ self.NU*3/(2*self.DX**3)*(
                           u[ 3:  , :, : ] - 2*u[ 2:-1, :, : ]
                       + 2*u[ 1:-2, :, : ] -   u[  :-3, :, : ]
                                          )
)
        return src
    # end calc_source_all_terms_3d

    def relax_pressure_poisson_3d(self, p,p_old, src):
        """Solves the poisson equation for a 3D pressure field
        by central differencing in both dimensions. This solves
        the laplace equation for a 3D pressure field when src=0
        NOTE: if p is initialized to all zeros, then the first
        calculation of pct_err will raise a divide by zero Runtime
        Warning; however, since inf > self.MAE that is okay and it
        functions correctly."""
            
            ##now, in x
        p[ 1:-1, 1:-1, 1:-1 ] = (self.DY**2*self.DZ**2*(
            p_old[ 2:  , 1:-1, 1:-1 ] + p_old[  :-2, 1:-1, 1:-1 ])
            / (2*(self.DY**2*self.DZ**2 + self.DX**2*self.DZ**2
                                        + self.DX**2*self.DY**2)))
            ##next, in y
        p[ 1:-1, 1:-1, 1:-1 ] += (self.DX**2*self.DZ**2*(
            p_old[ 1:-1, 2:  , 1:-1 ] + p_old[ 1:-1,  :-2, 1:-1 ])
            / (2*(self.DY**2*self.DZ**2 + self.DX**2*self.DZ**2
                                        + self.DX**2*self.DY**2)))
            ##last, in z
        p[ 1:-1, 1:-1, 1:-1 ] += (self.DX**2*self.DY**2*(
            p_old[ 1:-1, 1:-1, 2:   ] + p_old[ 1:-1, 1:-1,  :-2 ])
            / (2*(self.DY**2*self.DZ**2 + self.DX**2*self.DZ**2
                                        + self.DX**2*self.DY**2)))

            ##last, the source
        p[ 1:-1, 1:-1, 1:-1 ]  -= (self.DX**2*self.DY**2*self.DZ**2
                                     * src[ 1:-1, 1:-1, 1:-1 ]
            / (2*(self.DY**2*self.DZ**2 + self.DX**2*self.DZ**2
                                        + self.DX**2*self.DY**2)))
        return p
    # end relax_pressure_poisson_3d
    
# end class The_Fluid_Solver_3D

## end the_fluid_solver_3d.py
