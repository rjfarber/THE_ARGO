## the_bnds_setter_2d.py
## by Ryan Farber 20 January  2015
## Last modified: 14 February 2015

from the_solver_2d import The_Solver_2D

class The_Bnds_Setter_2D(The_Solver_2D):
    """
    This file contains all the boundary condition functions to be used by the_argo. Currently available boundary conditions are

    - Fixed       (Dirichlet or Neumann depending upon input)
    - Periodic

    and reflecting boundaries will hopefully one day be implemented.

    For Neumann Fixed boundary conditions, set the walls BEFORE the corners; since the walls overwrite the corner values.
    """
    def set_bnds_fixed_2d(self, f, amt):
        """Sets all walls of input 2d numpy array to a fixed amount of type int or float for Dirichlet boundary conditions."""

        f[ 0, : ] = amt; f[ -1,  : ] = amt
        f[ :, 0 ] = amt; f[  :, -1 ] = amt

        return f
    # end set_bnds_fixed_2d
    
    def set_bnds_fixed_2dXI(self, f, amt):
        """Sets initial wall along x of input 2d numpy array to a fixed amount of type int or float for Dirichlet (value) or Neumann (zero flux) or periodic boundary conditions."""

        f[ 0, : ] = amt

        return f
    # end set_bnds_fixed_2dXI

    def set_bnds_fixed_2dXF(self, f, amt):
        """Sets final wall along x of input 2d numpy array to a fixed amount of type int or float for Dirichlet (value) or Neumann (zero flux) boundary conditions."""

        f[ -1,  : ] = amt

        return f
    # end set_bnds_fixed_2dXF

    def set_bnds_fixed_2dYI(self, f, amt):
        """Sets initial wall along y of input 2d numpy array to a fixed amount of type int or float for Dirichlet (value) or Neumann (zero flux) or periodic boundary conditions."""

        f[ :, 0 ] = amt

        return f
    # end set_bnds_fixed_2dYI

    def set_bnds_fixed_2dYF(self, f, amt):
        """Sets final wall along y of input 2d numpy array to a fixed amount of type int or float for Dirichlet (value) or Neumann (zero flux) boundary conditions."""

        f[ :, -1 ] = amt

        return f
    # end set_bnds_fixed_2dYF

    def set_bnds_zero_flux_walls_2dX(self, f):
        """Sets zero flux boundary conditions on the walls of a 2D field."""

        f[ 0, : ] = -f[ 1, : ];  f[ -1,  : ] = -f[ -2,  : ]
        f[ :, 0 ] =  f[ :, 1 ];  f[  :, -1 ] =  f[  :, -2 ]

        return f
    # end set_bnds_zero_flux_walls_2dX

    def set_bnds_zero_flux_walls_2dY(self, f):
        """Sets zero flux boundary conditions on the walls of a 2D field."""

        f[ 0, : ] =  f[ 1, : ];  f[ -1,  : ] =  f[ -2,  : ]
        f[ :, 0 ] = -f[ :, 1 ];  f[  :, -1 ] = -f[  :, -2 ]

        return f
    # end set_bnds_zero_flux_walls_2dY

    def set_bnds_zero_flux_walls_scalar_2d(self, f):
        """Sets zero flux boundary conditions on the walls of a 2D field."""

        f[ 0, : ] =  f[ 1, : ];  f[ -1,  : ] =  f[ -2,  : ]
        f[ :, 0 ] =  f[ :, 1 ];  f[  :, -1 ] =  f[  :, -2 ]

        return f
    # end set_bnds_zero_flux_walls_scalar_2d

    def set_bnds_zero_flux_corners_2d(self, f):
        """Sets zero flux boundary conditions on the corners of a 2D field."""
        
        f[  0,  0 ] = (f[  1,  0 ] + f[  0,  1 ]) / 2.0    
        f[ -1,  0 ] = (f[ -2,  0 ] + f[ -1,  1 ]) / 2.0
        f[  0, -1 ] = (f[  1, -1 ] + f[  0, -2 ]) / 2.0
        f[ -1, -1 ] = (f[ -2, -1 ] + f[ -1, -2 ]) / 2.0
        
        return f
    # end set_bnds_zero_flux_corners_2d

    def set_bnds_periodic_diffusion_2dX(self, g):
        """Sets periodic boundary conditions along the x-dimension for any component of a 2D field (f) specifically for explicit central differencing diffusion of the input field (g) which has values of f at the old time step."""

        f  = self.NU*(self.DT/self.DX**2
    * (g[  0, 1:-1 ] - 2*g[ -1, 1:-1 ] + g[ -2, 1:-1 ]))
                                        
        f += self.NU*(self.DT/self.DY**2
    * (g[ -1, 2:   ] - 2*g[ -1, 1:-1 ] + g[ -1,  :-2 ]))
        
        return f
    # end set_bnds_periodic_diffusion_2dX

    def set_bnds_periodic_diffusion_2dY(self, g):
        """Sets periodic boundary conditions along the y-dimension for any component of a 2D field (f) specifically for explicit central differencing diffusion of the input field (g) which has values of f at the old time step."""

        f  = self.NU*(self.DT/self.DX**2
    * (g[ 2:  , -1 ] - 2*g[ 1:-1, -1 ] + g[  :-2, -1 ]))
            
        f += self.NU*(self.DT/self.DY**2
    * (g[ 1:-1,  0 ] - 2*g[ 1:-1, -1 ] + g[ 1:-1, -2 ]))
                              
        return f
    # end set_bnds_periodic_diffusion_2dY

    def set_bnds_periodic_diffusion_corners_2d(self, f,g):
        """Sets periodic boundary conditions of the corners for any component of a 2D field (f) specifically for explicit central differencing diffusion of the input field (g) which has values of f at the old time step."""
        f[  0,  0 ] += self.DT * (self.NU *
       (g[  1,  0 ] - 2*g[  0,  0 ] + g[ -1,  0 ]) / self.DX**2
                                 )
        f[  0,  0 ] += self.DT * (self.NU *
       (g[  0,  1 ] - 2*g[  0,  0 ] + g[  0, -1 ]) / self.DY**2
                                 )
        
        f[ -1,  0 ] += self.DT * (self.NU *
       (g[  0,  0 ] - 2*g[ -1,  0 ] + g[ -2,  0 ]) / self.DX**2
                                 )                                   
        f[ -1,  0 ] += self.DT * (self.NU *
       (g[ -1,  1 ] - 2*g[ -1,  0 ] + g[ -1, -1 ]) / self.DY**2
                                 )
        
        f[  0, -1 ] += self.DT * (self.NU *
       (g[  1, -1 ] - 2*g[  0, -1 ] + g[ -1, -1 ]) / self.DX**2
                                 )
        f[  0, -1 ] += self.DT * (self.NU *
       (g[  0,  0 ] - 2*g[  0, -1 ] + g[  0, -2 ]) / self.DY**2
                                 )

        f[ -1, -1 ] += self.DT * (self.NU *
       (g[  0, -1 ] - 2*g[ -1, -1 ] + g[ -2, -1 ]) / self.DX**2
                                 )                                   
        f[ -1, -1 ] += self.DT * (self.NU *
       (g[ -1,  0 ] - 2*g[ -1, -1 ] + g[ -1, -2 ]) / self.DY**2
                                 )
        return f
    # end set_bnds_periodic_diffusion_corners_2d

    def set_bnds_periodic_apply_pressure_2dX(self, fx,fy, p):
        """Sets periodic boundary conditions along the x-dimension for a 2D field (f) specifically for a pressure gradient."""

        fx[ -1 , 1:-1 ] -= self.DT/(2*self.RHO*self.DX)*(
                            p[  0 , 1:-1 ] - p[  -2 , 1:-1 ])

        fy[ -1 , 1:-1 ] -= self.DT/(2*self.RHO*self.DY)*(
                            p[ -1 , 2:   ] - p[  -1 ,  :-2 ])
        return [fx, fy]
    # end set_bnds_periodic_apply_pressure_2dX

    def set_bnds_periodic_apply_pressure_2dY(self, fx,fy, p):
        """Sets periodic boundary conditions along the y-dimension for a 2D field (f) specifically for a pressure gradient."""

        fx[ 1:-1,  -1 ] -= self.DT/(2*self.RHO*self.DX)*(
                            p[ 2:  ,  -1 ] - p[  :-2,  -1 ])

        fy[ 1:-1,  -1 ] -= self.DT/(2*self.RHO*self.DY)*(
                            p[ 1:-1,   0 ] - p[ 1:-1 , -2 ])
        return [fx, fy]
    # end set_bnds_periodic_apply_pressure_2dY

    def set_bnds_periodic_apply_force_2dX(self, f,g):
        """Sets periodic boundary conditions along the x-dimension for a 2D field (f) specifically for the acceleration due to an external force."""

        f[ -1 , 1:-1 ] += self.DT * g[ -1, 1:-1 ]

        return f
    # end set_bnds_periodic_apply_force_2dX

    def set_bnds_periodic_apply_force_2dY(self, f,g):
        """Sets periodic boundary conditions along the y-dimension for a 2D field (f) specifically for the acceleration due to an external force."""
            
        f[ 1:-1, -1 ] += self.DT * g[ 1:-1, -1 ]

        return f
    # end set_bnds_periodic_apply_force_2dY

    def set_bnds_periodic_src_2dX(self, src, u,v):
        """Sets boundary conditions along the x-dimension for the source term of the poisson eqn of velocity in the navier-stokes equation."""
        
        src[ -1, 1:-1 ] = self.RHO*(1.0/self.DT*(
            (u[  0, 1:-1 ] - u[ -2, 1:-1 ])/(2*self.DX)
      +     (v[ -1, 2:   ] - v[ -1,  :-2 ])/(2*self.DY)
                                      )
      -    ((u[  0, 1:-1 ] - u[ -2, 1:-1 ]) / (2*self.DX))**2
      -  2*((u[ -1, 2:   ] - u[ -1,  :-2 ]) / (2*self.DY)
      *     (v[  0, 1:-1 ] - v[ -2, 1:-1 ]) / (2*self.DX))
      -    ((v[ -1, 2:   ] - v[ -1,  :-2 ]) / (2*self.DY))**2
                                   )
        ##set first equal to last
        src[ 0, 1:-1 ] = src[ -1, 1:-1 ]

        return src
    # end set_bnds_periodic_src_2dX

    def set_bnds_periodic_src_2dY(self, src, u,v):
        """Sets boundary conditions along the y-dimension for the source term of the poisson eqn of velocity in the navier-stokes equation."""
        
        src[ 1:-1, -1 ] = self.RHO*(1.0/self.DT*(
            (u[ 2:  , -1 ] - u[  :-2, -1 ]) / (2*self.DX)
      +     (v[ 1:-1,  0 ] - v[ 1:-1, -2 ]) / (2*self.DY)
                                                )
      -    ((u[ 2:  , -1 ] - u[  :-2, -1 ]) / (2*self.DX))**2
      -  2*((u[ 1:-1,  0 ] - u[ 1:-1, -2 ]) / (2*self.DY)
      *     (v[ 2:  , -1 ] - v[  :-2, -1 ]) / (2*self.DX))
      -    ((v[ 1:-1,  0 ] - v[ 1:-1, -2 ]) / (2*self.DY))**2
                                   )
        ##set first equal to last
        src[ 1:-1, 0 ] = src[ 1:-1, -1 ]

        return src
    # end set_bnds_periodic_src_2dY

    def set_bnds_periodic_src_div_only_2dX(self, src, u,v):
        """Sets boundary conditions along the x-dimension for the source term of the poisson eqn of velocity in the navier-stokes equation."""
        
        src[ -1, 1:-1 ] = self.RHO*(1.0/self.DT*(
            (u[  0, 1:-1 ] - u[ -2, 1:-1 ])/(2*self.DX)
      +     (v[ -1, 2:   ] - v[ -1,  :-2 ])/(2*self.DY)
                                      )
                               )
        ##set first equal to last
        src[ 0, 1:-1 ] = src[ -1, 1:-1 ]

        return src
    # end set_bnds_periodic_src_div_only_2dX

    def set_bnds_periodic_src_div_only_2dY(self, src, u,v):
        """Sets boundary conditions along the y-dimension for the source term of the poisson eqn of velocity in the navier-stokes equation."""
        
        src[ 1:-1, -1 ] = self.RHO*(1.0/self.DT*(
            (u[ 2:  , -1 ] - u[  :-2, -1 ]) / (2*self.DX)
      +     (v[ 1:-1,  0 ] - v[ 1:-1, -2 ]) / (2*self.DY)
                                                )
                               )
        ##set first equal to last
        src[ 1:-1, 0 ] = src[ 1:-1, -1 ]

        return src
    # end set_bnds_periodic_src_div_only_2dY

    def set_bnds_periodic_pressure_2dX(self, p,p_old, src):
        """Sets periodic boundary conditions along the x-dimension for pressure for the navier-stokes equations."""
            ## first in x
        p[ -1, 1:-1 ]  = ( self.DY**2 /
( 2*(self.DX**2 + self.DY**2) ) *
(p_old[  0, 1:-1 ] + p_old[ -2, 1:-1 ]) )
                           
            ## now in y
        p[ -1, 1:-1 ] += ( self.DX**2 /
( 2*(self.DX**2 + self.DY**2) )
* (p_old[ -1, 2:   ] + p_old[ -1,  :-2 ]) )
                           
            ## last, the source
        p[ -1, 1:-1 ] -= (src[ -1, 1:-1 ] *
self.DX**2 * self.DY**2
/ (2*(self.DX**2 + self.DY**2)) )
        
            ##set first equal to last
        p[ 0, 1:-1 ] = p[ -1, 1:-1 ]

        return p
    # end set_bnds_periodic_pressure_2dX

    def set_bnds_periodic_pressure_2dY(self, p,p_old, src):
        """Sets periodic boundary conditions along the y-dimension for
        pressure for the navier-stokes equations."""
            ## first in x
        p[ 1:-1, -1 ]  = (self.DY**2 / ( 2*(self.DX**2 + self.DY**2) )
                       * (p_old[ 2:  , -1 ] + p_old[  :-2, -1 ]) )
            ## now in y
        p[ 1:-1, -1 ] += (self.DX**2 / ( 2*(self.DX**2 + self.DY**2)) 
                       * (p_old[ 1:-1,  0 ] + p_old[ 1:-1, -2 ]) )
                           
            ## last, the source
        p[ 1:-1, -1 ] -= (src[ 1:-1, -1 ] * self.DX**2 * self.DY**2
                           / (2*(self.DX**2 + self.DY**2)) )
        
            ##set first equal to last
        p[ 1:-1, 0 ] = p[ 1:-1, -1 ]

        return p
    # end set_bnds_periodic_pressure_2dY

# end The_Bnds_Setter_2D

## end the_bnds_setter_2d.py
