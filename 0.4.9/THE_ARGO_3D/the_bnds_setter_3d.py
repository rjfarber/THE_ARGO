## the_bnds_setter_3d.py
## by Ryan Farber 20 January 2015
## Last modified: 27 January 2015
"""
This file contains all the boundary condition functions to be used by
the_argo. Currently available boundary conditions are

- Zero Flux   (Neumann)
- Fixed       (Dirichlet)
- Periodic

and reflecting boundaries will hopefully one day be implemented.
"""

##3D ZERO FLUX BOUNDARY CONDITIONS
def set_bnds_zero_flux_3dX(fx):
    """Sets zero flux boundary conditions for the x-component
    of a 3D (velocity) field."""
    fx[0,:,:] = -fx[1,:,:];  fx[ -1,  :,  : ] = -fx[ -2,  :,  : ]
    fx[:,0,:] =  fx[:,1,:];  fx[  :, -1,  : ] =  fx[  :, -2,  : ]
    fx[:,:,0] =  fx[:,:,1];  fx[  :,  :, -1 ] =  fx[  :,  :, -2 ]

    return fx
# end set_bnds_zero_flux_3dX

def set_bnds_zero_flux_3dY(fy):
    """Sets zero flux boundary conditions for the y-component
    of a 3D (velocity) field."""

    fy[0,:,:] =  fy[1,:,:];  fy[ -1,  :,  : ] =  fy[ -2,  :,  : ]
    fy[:,0,:] = -fy[:,1,:];  fy[  :, -1,  : ] = -fy[  :, -2,  : ]
    fy[:,:,0] =  fy[:,:,1];  fy[  :,  :, -1 ] =  fy[  :,  :, -2 ]

    return fy
# end set_bnds_zero_flux_3dY

def set_bnds_zero_flux_3dZ(fz):
    """Sets zero flux boundary conditions for the z-component
    of a 3D (velocity) field."""

    fz[0,:,:] =  fz[1,:,:];  fz[ -1,  :,  : ] =  fz[ -2,  :,  : ]
    fz[:,0,:] =  fz[:,1,:];  fz[  :, -1,  : ] =  fz[  :, -2,  : ]
    fz[:,:,0] = -fz[:,:,1];  fz[  :,  :, -1 ] = -fz[  :,  :, -2 ]

    return fz
# end set_bnds_zero_flux_3dZ

def set_bnds_zero_flux_scalar_3d(f):
    """Sets zero flux boundary conditions for a scalar
    3D (pressure) field."""

    f[0,:,:] = f[1,:,:];  f[ -1,  :,  : ] = f[ -2,  :,  : ]
    f[:,0,:] = f[:,1,:];  f[  :, -1,  : ] = f[  :, -2,  : ]
    f[:,:,0] = f[:,:,1];  f[  :,  :, -1 ] = f[  :,  :, -2 ]

    return f
# end set_bnds_zero_flux_scalar_3d

def set_bnds_zero_flux_corners_3d(f):
    """Sets zero flux boundary conditions on the corners of
    a 3D field."""
    f[  0,  0,  0 ] = (f[  1,  0,  0 ] + f[  0,  1,  0 ]
                             + f[  0,  0,  1 ]) / 3.0
    
    f[ -1,  0,  0 ] = (f[ -2,  0,  0 ] + f[ -1,  1,  0 ]
                             + f[ -1,  0,  1 ]) / 3.0
    
    f[  0, -1,  0 ] = (f[  1, -1,  0 ] + f[  0, -2,  0 ]
                             + f[  0, -1,  1 ]) / 3.0

    f[  0,  0, -1 ] = (f[  1,  0, -1 ] + f[  0,  1, -1 ]
                             + f[  0,  0, -2 ]) / 3.0

    f[ -1, -1,  0 ] = (f[ -2, -1,  0 ] + f[ -1, -2,  0 ]
                             + f[ -1, -1,  1 ]) / 3.0

    f[  0, -1, -1 ] = (f[  1, -1, -1 ] + f[  0, -2, -1 ]
                             + f[  0, -1, -2 ]) / 3.0

    f[ -1,  0, -1 ] = (f[ -2,  0, -1 ] + f[ -1,  1, -1 ]
                             + f[ -1,  0, -2 ]) / 3.0

    f[ -1, -1, -1 ] = (f[ -2, -1, -1 ] + f[ -1, -2, -1 ]
                             + f[ -1, -1, -2 ]) / 3.0
    return f
# end set_bnds_zero_flux_corners_3d




"*********************************************************************"
#*********************************************************************#




##3D FIXED BOUNDARY CONDITIONS
def set_bnds_fixed_3dX(f, amtI, amtF):
    """Sets faces to a fixed (input) value 'amt' along the
    x-dimension; Dirichlet boundary conditions."""

    f[ 0, :, : ] = amtI; f[ -1,  :,  : ] = amtF

    return f
# end set_bnds_fixed_3dX

def set_bnds_fixed_3dY(f, amtI, amtF):
    """Sets faces to a fixed (input) value 'amt' along the
    y-dimension; Dirichlet boundary conditions."""

    f[ :, 0, : ] = amtI; f[  :, -1,  : ] = amtF

    return f
# end set_bnds_fixed_3dY

def set_bnds_fixed_3dZ(f, amtI, amtF):
    """Sets faces to a fixed (input) value 'amt' along the
    z-dimension; Dirichlet boundary conditions."""

    f[ :, :, 0 ] = amtI; f[  :,  :, -1 ] = amtF

    return f
# end set_bnds_fixed_3dZ




"*********************************************************************"
#*********************************************************************#




##3D PERIODIC BOUNDARY CONDITIONS
def set_bnds_periodic_3dX(f):
    """Completes periodic boundary conditions along the x-dimension.
    Sets first equal to the last."""
    
    f[ 0, :, : ] = f[ -1,  :,  : ]

    return f
# end set_bnds_periodic_3dX

def set_bnds_periodic_3dY(f):
    """Completes periodic boundary conditions along the y-dimension.
    Sets first equal to the last."""
    
    f[ :, 0, : ] = f[  :, -1,  : ]

    return f
# end set_bnds_periodic_3dY

def set_bnds_periodic_3dZ(f):
    """Completes periodic boundary conditions along the z-dimension.
    Sets first equal to the last."""
    
    f[ :, :, 0 ] = f[  :,  :, -1 ]

    return f
# end set_bnds_periodic_3dZ

def set_bnds_periodic_diffusion_3dX(DT,DX,DY,DZ,NU, f,f_old):
    """Sets periodic boundary conditions along the x-dimension for
    any component of a 3D field (f) specifically for diffusion."""

    f[ -1, 1:-1, 1:-1 ] += NU*((DT/DX**2)
* (f_old[  0, 1:-1, 1:-1 ] - 2*f_old[ -1, 1:-1, 1:-1 ]
+  f_old[ -2, 1:-1, 1:-1 ]))
                                    
    f[ -1, 1:-1, 1:-1 ] += NU*((DT/DY**2)
* (f_old[ -1, 2:  , 1:-1 ] - 2*f_old[ -1, 1:-1, 1:-1 ]
+  f_old[ -1,  :-2, 1:-1 ]))

    f[ -1, 1:-1, 1:-1 ] += NU*((DT/DZ**2)
* (f_old[ -1, 1:-1, 2:   ] - 2*f_old[ -1, 1:-1, 1:-1 ]
+  f_old[ -1, 1:-1,  :-2 ]))

    return f
# end set_bnds_periodic_diffusion_3dX

def set_bnds_periodic_diffusion_3dY(DT,DX,DY,DZ,NU, f,f_old):
    """Sets periodic boundary conditions along the y-dimension for
    any component of a 3D field (f) specifically for diffusion."""

    f[ 1:-1, -1, 1:-1 ] += NU*((DT/DX**2)
* (f_old[ 2:  , -1, 1:-1 ] - 2*f_old[ 1:-1, -1, 1:-1 ]
+  f_old[  :-2, -1, 1:-1 ]))
        
    f[ 1:-1, -1, 1:-1 ] += NU*((DT/DY**2)
* (f_old[ 1:-1,  0, 1:-1 ] - 2*f_old[ 1:-1, -1, 1:-1 ]
+  f_old[ 1:-1, -2, 1:-1 ]))

    f[ 1:-1, -1, 1:-1 ] += NU*((DT/DZ**2)
* (f_old[ 1:-1, -1, 2:   ] - 2*f_old[ 1:-1, -1, 1:-1 ]
+  f_old[ 1:-1, -1,  :-2 ]))

    return f
# end set_bnds_periodic_diffusion_3dY

def set_bnds_periodic_diffusion_3dZ(DT,DX,DY,DZ,NU, f,f_old):
    """Sets periodic boundary conditions along the z-dimension for
    any component of a 3D field (f) specifically for diffusion."""

    f[ 1:-1, 1:-1, -1 ] += NU*((DT/DX**2)
* (f_old[ 2:  , 1:-1, -1 ] - 2*f_old[ 1:-1, 1:-1, -1 ]
+  f_old[  :-2, 1:-1, -1 ]))
        
    f[ 1:-1, 1:-1, -1 ] += NU*((DT/DY**2)
* (f_old[ 1:-1, 2:  , -1 ] - 2*f_old[ 1:-1, 1:-1, -1 ]
+  f_old[ 1:-1,  :-2, -1 ]))

    f[ 1:-1, 1:-1, -1 ] += NU*((DT/DZ**2)
* (f_old[ 1:-1, 1:-1,  0 ] - 2*f_old[ 1:-1, 1:-1, -1 ]
+  f_old[ 1:-1, 1:-1, -2 ]))

    return f
# end set_bnds_periodic_diffusion_3dZ

def set_bnds_periodic_apply_pressure_3dX(DT,DX,DY,DZ,RHO, fx,fy,fz, p):
    """Sets periodic boundary conditions along the x-dimension for
    a 3D field (f) specifically for a pressure gradient."""

    fx[  -1 , 1:-1, 1:-1 ] -= DT/(2*RHO*DX)*(
        p[   0 , 1:-1, 1:-1 ] - p[  -2 , 1:-1, 1:-1 ])

    fy[  -1 , 1:-1, 1:-1 ] -= DT/(2*RHO*DY)*(
        p[  -1 , 2:  , 1:-1 ] - p[  -1 ,  :-2, 1:-1 ])

    fz[  -1 , 1:-1, 1:-1 ] -= DT/(2*RHO*DZ)*(
        p[  -1 , 1:-1, 2:   ] - p[  -1 , 1:-1,  :-2 ])

    return [fx, fy, fz]
# end set_bnds_periodic_apply_pressure_3dX

def set_bnds_periodic_apply_pressure_3dY(DT,DX,DY,DZ,RHO, fx,fy,fz, p):
    """Sets periodic boundary conditions along the y-dimension for
    a 3D field (f) specifically for a pressure gradient."""

    fx[ 1:-1,  -1 , 1:-1 ] -= DT/(2*RHO*DX)*(
        p[ 2:  ,  -1 , 1:-1 ] - p[  :-2,  -1 , 1:-1 ])

    fy[ 1:-1,  -1 , 1:-1 ] -= DT/(2*RHO*DY)*(
        p[ 1:-1,   0 , 1:-1 ] - p[  -1 ,  -2 , 1:-1 ])

    fz[ 1:-1,  -1 , 1:-1 ] -= DT/(2*RHO*DZ)*(
        p[ 1:-1,  -1 , 2:   ] - p[ 1:-1,  -1 ,  :-2 ])

    return [fx, fy, fz]
# end set_bnds_periodic_apply_pressure_3dY

def set_bnds_periodic_apply_pressure_3dZ(DT,DX,DY,DZ,RHO, fx,fy,fz, p):
    """Sets periodic boundary conditions along the z-dimension for
    a 3D field (f) specifically for a pressure gradient."""

    fx[ 1:-1, 1:-1,  -1  ] -= DT/(2*RHO*DX)*(
        p[ 2:  , 1:-1,  -1  ] - p[  :-2, 1:-1,  -1  ])

    fy[ 1:-1, 1:-1,  -1  ] -= DT/(2*RHO*DY)*(
        p[ 1:-1, 2:  ,  -1  ] - p[ 1:-1,  :-2,  -1 ])

    fz[ 1:-1, 1:-1,  -1  ] -= DT/(2*RHO*DZ)*(
        p[ 1:-1, 1:-1,   0  ] - p[ 1:-1, 1:-1,  -2 ])

    return [fx, fy, fz]
# end set_bnds_periodic_apply_pressure_3dZ

def set_bnds_periodic_apply_force_3dX(DT, f,g):
    """Sets periodic boundary conditions along the x-dimension
    for a 3D field (f) specifically for the acceleration due to
    an external force."""
        ##along x
    f[  -1 , 1:-1, 1:-1 ] += DT * g[  -1 , 1:-1, 1:-1 ]

    return f
# end set_bnds_periodic_apply_force_3dX

def set_bnds_periodic_apply_force_3dY(DT, f,g):
    """Sets periodic boundary conditions along the y-dimension
    for a 3D field (f) specifically for the acceleration due to
    an external force."""
        ##along y
    f[ 1:-1,  -1 , 1:-1 ] += DT * g[ 1:-1,  -1 , 1:-1 ]

    return f
# end set_bnds_periodic_apply_force_3dY

def set_bnds_periodic_apply_force_3dZ(DT, f,g):
    """Sets periodic boundary conditions along the z-dimension
    for a 3D field (f) specifically for the acceleration due to
    an external force."""
        ##along z
    f[ 1:-1, 1:-1,  -1  ] += DT * g[ 1:-1, 1:-1,  -1  ]

    return f
# end set_bnds_periodic_apply_force_3dZ

##3D PRESSURE BOUNDARY CONDITIONS    
def set_bnds_pressure_3d(p, tag):
    """Sets boundary conditions for a specific problem."""

    if   tag == 'lid' or tag[:6] == 'cavity':
            ##dp/dy = 0 @ y = 0,2; dp/dx = 0 @ x = 0; p = 0 @ x = 2
            ##(in the order listed above)
        p[ -1,  :, : ] = p[ -2, :, : ]
        p[  0,  :, : ] = p[  1, :, : ]
        p[  :,  0, : ] = p[  :, 1, : ]
        p[  :, -1, : ] = 0.0
        
    else:
            ##Default: Fixed at faces
        p[ 0, :, : ] = 0.0; p[ -1,  :,  : ] = 0.0
        p[ :, 0, : ] = 0.0; p[  :, -1,  : ] = 0.0
        p[ :, :, 0 ] = 0.0; p[  :,  :, -1 ] = 0.0

    return p
# end set_bnds_pressure_3d

def set_bnds_periodic_pressure_3dX(DX,DY,DZ, p,p_old, src):
    """Sets periodic boundary conditions along the x-dimension
    for pressure of a poisson eqn."""

        ##now, in x
    p[ -1, 1:-1, 1:-1 ]  = (DY**2*DZ**2*(
        p_old[ 0  , 1:-1, 1:-1 ] + p_old[  -2, 1:-1, 1:-1 ])
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))
        ##next, in y
    p[ -1, 1:-1, 1:-1 ] += (DX**2*DZ**2*(
        p_old[ -1, 2:  , 1:-1 ] + p_old[ -1,  :-2, 1:-1 ])
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))
        ##last, in z
    p[ -1, 1:-1, 1:-1 ] += (DX**2*DY**2*(
        p_old[ -1, 1:-1, 2:   ] + p_old[ -1, 1:-1,  :-2 ])
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))

        ##last, the source
    p[ -1, 1:-1, 1:-1 ]  -= (DX**2*DY**2*DZ**2 * src[ -1, 1:-1, 1:-1 ]
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))
        ##set first equal to last
    p[ 0, :, : ] = p[ -1, :, : ]
    
    return p
# end set_bnds_periodic_pressure_3dX

def set_bnds_periodic_pressure_3dY(DX,DY,DZ, p,p_old, src):
    """Sets periodic boundary conditions along the y-dimension
    for pressure of a poisson eqn."""

        ##now, in x
    p[ 1:-1, -1, 1:-1 ]  = (DY**2*DZ**2*(
        p_old[ 2:  , -1, 1:-1 ] + p_old[  :-2, -1, 1:-1 ])
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))
        ##next, in y
    p[ 1:-1, -1, 1:-1 ] += (DX**2*DZ**2*(
        p_old[ 1:-1,  0, 1:-1 ] + p_old[ 1:-1, -2, 1:-1 ])
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))
        ##last, in z
    p[ 1:-1, -1, 1:-1 ] += (DX**2*DY**2*(
        p_old[ 1:-1, -1, 2:   ] + p_old[ 1:-1, -1,  :-2 ])
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))

        ##last, the source
    p[ 1:-1, -1, 1:-1 ]  -= (DX**2*DY**2*DZ**2 * src[ 1:-1, -1, 1:-1 ]
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))
        ##set first equal to last
    p[ :, 0, : ] = p[ :, -1, : ]
    
    return p
# end set_bnds_periodic_pressure_3dY

def set_bnds_periodic_pressure_3dZ(DX,DY,DZ, p,p_old, src):
    """Sets periodic boundary conditions along the z-dimension
    for pressure of a poisson eqn."""

        ##now, in x
    p[ 1:-1, 1:-1, -1 ]  = (DY**2*DZ**2*(
        p_old[ 2:  , 1:-1, -1 ] + p_old[  :-2, 1:-1, -1 ])
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))
        ##next, in y
    p[ 1:-1, 1:-1, -1 ] += (DX**2*DZ**2*(
        p_old[ 1:-1, 1:-1, -1 ] + p_old[ 1:-1, 1:-1, -1 ])
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))
        ##last, in z
    p[ 1:-1, 1:-1, -1 ] += (DX**2*DY**2*(
        p_old[ 1:-1, 1:-1,  0 ] + p_old[ 1:-1, 1:-1, -2 ])
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))

        ##last, the source
    p[ 1:-1, 1:-1, -1 ]  -= (DX**2*DY**2*DZ**2 * src[ 1:-1, 1:-1, -1 ]
            / (2*(DY**2*DZ**2 + DX**2*DZ**2 + DX**2*DY**2)))
        ##set first equal to last
    p[ :, :, 0 ] = p[ :, :, -1 ]
    
    return p
# end set_bnds_periodic_pressure_3dZ

def set_bnds_periodic_src_3dX(DT,DX,DY,DZ, RHO, src, u,v,w):
    """Sets boundary conditions along the x-dimension for the
    source term of the poisson eqn of velocity in the
    navier-stokes equation."""

    src[ -1, 1:-1, 1:-1 ] = RHO*(
1.0/DT*(
                (u[  0 , 1:-1, 1:-1 ] - u[ -2, 1:-1, 1:-1 ])
                                        / (2*DX)  +
                (v[ -1 , 2:  , 1:-1 ] - v[ -1,  :-2, 1:-1 ])
                                        / (2*DY)  +
                (w[ -1 , 1:-1, 2:   ] - w[ -1, 1:-1,  :-2 ])
                                        / (2*DZ)
        )
          -
            (   (u[  0 , 1:-1, 1:-1 ] - u[ -2, 1:-1, 1:-1 ])
                                        / (2*DX)
            )**2
          -
            (   (v[ -1 , 2:  , 1:-1 ] - v[ -1,  :-2, 1:-1 ])
                                        / (2*DY)
            )**2
          -
            (   (w[ -1 , 1:-1, 2:   ] - w[ -1, 1:-1,  :-2 ])
                                        / (2*DZ)
            )**2
          -
          2*(   (v[  0 , 1:-1, 1:-1 ] - v[ -2, 1:-1, 1:-1 ])
                                        / (2*DX)
              * (u[ -1 , 2:  , 1:-1 ] - u[ -1,  :-2, 1:-1 ])
                                        / (2*DY)
            )
          -
          2*(   (w[  0 , 1:-1, 1:-1 ] - w[ -2, 1:-1, 1:-1 ])
                                        / (2*DX)
              * (u[ -1 , 1:-1, 2:   ] - u[ -1, 1:-1,  :-2 ])
                                        / (2*DZ)
            )
          -
          2*(   (w[ -1 , 2:  , 1:-1 ] - w[ -1,  :-2, 1:-1 ])
                                        / (2*DY)
              * (v[ -1 , 1:-1, 2:   ] - v[ -1, 1:-1,  :-2 ])
                                        / (2*DZ)
            )
                                      )
        ##set first equal to last
    src[ 0, :, : ] = src[ -1, :, : ]
    
    return src
# end set_bnds_periodic_src_3dX

def set_bnds_periodic_src_3dY(DT,DX,DY,DZ, RHO, src, u,v,w):
    """Sets boundary conditions along the y-dimension for the
    source term of the poisson eqn of velocity in the
    navier-stokes equation."""

    src[ 1:-1, -1, 1:-1 ] = RHO*(
1.0/DT*(
                (u[ 2:  , -1, 1:-1 ] - u[  :-2, -1, 1:-1 ])
                                        / (2*DX)  +
                (v[ 1:-1,  0, 1:-1 ] - v[ 1:-1, -2, 1:-1 ])
                                        / (2*DY)  +
                (w[ 1:-1, -1, 2:   ] - w[ 1:-1, -1,  :-2 ])
                                        / (2*DZ)
        )
          -
            (   (u[ 2:  , -1, 1:-1 ] - u[  :-2, -1, 1:-1 ])
                                        / (2*DX)
            )**2
          -
            (   (v[ 1:-1,  0, 1:-1 ] - v[ 1:-1, -2, 1:-1 ])
                                        / (2*DY)
            )**2
          -
            (   (w[ 1:-1, -1, 2:   ] - w[ 1:-1, -1,  :-2 ])
                                        / (2*DZ)
            )**2
          -
          2*(   (v[ 2:  , -1, 1:-1 ] - v[  :-2, -1, 1:-1 ])
                                        / (2*DX)
              * (u[ 1:-1,  0, 1:-1 ] - u[ 1:-1, -2, 1:-1 ])
                                        / (2*DY)
            )
          -
          2*(   (w[ 2:  , -1, 1:-1 ] - w[  :-2, -1, 1:-1 ])
                                        / (2*DX)
              * (u[ 1:-1, -1, 2:   ] - u[ 1:-1, -1,  :-2 ])
                                        / (2*DZ)
            )
          -
          2*(   (w[ 1:-1,  0, 1:-1 ] - w[ 1:-1, -2, 1:-1 ])
                                        / (2*DY)
              * (v[ 1:-1, -1, 2:   ] - v[ 1:-1, -1,  :-2 ])
                                        / (2*DZ)
            )
                                      )
        ##set first equal to last
    src[ :, 0, : ] = src[ :, -1, : ]
    
    return src
# end set_bnds_periodic_src_3dY

def set_bnds_periodic_src_3dZ(DT,DX,DY,DZ, RHO, src, u,v,w):
    """Sets boundary conditions along the z-dimension for the
    source term of the poisson eqn of velocity in the
    navier-stokes equation."""

    src[ 1:-1, 1:-1, -1 ] = RHO*(
1.0/DT*(
                (u[ 2:  , 1:-1, -1 ] - u[  :-2, 1:-1, -1 ])
                                        / (2*DX)  +
                (v[ 1:-1, 2:  , -1 ] - v[ 1:-1,  :-2, -1 ])
                                        / (2*DY)  +
                (w[ 1:-1, 1:-1,  0 ] - w[ 1:-1, 1:-1, -2 ])
                                        / (2*DZ)
        )
          -
            (   (u[ 2:  , 1:-1, -1 ] - u[  :-2, 1:-1, -1 ])
                                        / (2*DX)
            )**2
          -
            (   (v[ 1:-1, 2:  , -1 ] - v[ 1:-1,  :-2, -1 ])
                                        / (2*DY)
            )**2
          -
            (   (w[ 1:-1, 1:-1,  0 ] - w[ 1:-1, 1:-1, -2 ])
                                        / (2*DZ)
            )**2
          -
          2*(   (v[ 2:  , 1:-1, -1 ] - v[  :-2, 1:-1, -1 ])
                                        / (2*DX)
              * (u[ 1:-1, 2:  , -1 ] - u[ 1:-1,  :-2, -1 ])
                                        / (2*DY)
            )
          -
          2*(   (w[ 2:  , 1:-1, -1 ] - w[  :-2, 1:-1, -1 ])
                                        / (2*DX)
              * (u[ 1:-1, 1:-1,  0 ] - u[ 1:-1, 1:-1, -2 ])
                                        / (2*DZ)
            )
          -
          2*(   (w[ 1:-1, 2:  , -1 ] - w[ 1:-1,  :-2, -1 ])
                                        / (2*DY)
              * (v[ 1:-1, 1:-1,  0 ] - v[ 1:-1, 1:-1, -2 ])
                                        / (2*DZ)
            )
                                      )
        ##set first equal to last
    src[ :, :, 0 ] = src[ :, :, -1 ]
    
    return src
# end set_bnds_periodic_src_3dZ

def set_bnds_periodic_src_div_only_3dX(DT,DX,DY,DZ, RHO, src, u,v,w):
    """Sets boundary conditions along the x-dimension for the
    source term of the poisson eqn of velocity in the
    navier-stokes equation, including only the divergence term."""

    src[ -1, 1:-1, 1:-1 ] = RHO*(
1.0/DT*(
                (u[  0 , 1:-1, 1:-1 ] - u[ -2, 1:-1, 1:-1 ])
                                        / (2*DX)  +
                (v[ -1 , 2:  , 1:-1 ] - v[ -1,  :-2, 1:-1 ])
                                        / (2*DY)  +
                (w[ -1 , 1:-1, 2:   ] - w[ -1, 1:-1,  :-2 ])
                                        / (2*DZ)
       )
                                )
        ##set first equal to last
    src[ 0, :, : ] = src[ -1, :, : ]
    
    return src
# end set_bnds_periodic_src_div_only_3dX

def set_bnds_periodic_src_div_only_3dY(DT,DX,DY,DZ, RHO, src, u,v,w):
    """Sets boundary conditions along the y-dimension for the
    source term of the poisson eqn of velocity in the
    navier-stokes equation, including only the divergence term."""

    src[ 1:-1, -1, 1:-1 ] = RHO*(
1.0/DT*(
                (u[ 2:  , -1, 1:-1 ] - u[  :-2, -1, 1:-1 ])
                                        / (2*DX)  +
                (v[ 1:-1,  0, 1:-1 ] - v[ 1:-1, -2, 1:-1 ])
                                        / (2*DY)  +
                (w[ 1:-1, -1, 2:   ] - w[ 1:-1, -1,  :-2 ])
                                        / (2*DZ)
        )
                                )
        ##set first equal to last
    src[ :, 0, : ] = src[ :, -1, : ]
    
    return src
# end set_bnds_periodic_src_div_only_3dY

def set_bnds_periodic_src_div_only_3dZ(DT,DX,DY,DZ, RHO, src, u,v,w):
    """Sets boundary conditions along the z-dimension for the
    source term of the poisson eqn of velocity in the
    navier-stokes equation, including only the divergence term."""

    src[ 1:-1, 1:-1, -1 ] = RHO*(
1.0/DT*(
                (u[ 2:  , 1:-1, -1 ] - u[  :-2, 1:-1, -1 ])
                                        / (2*DX)  +
                (v[ 1:-1, 2:  , -1 ] - v[ 1:-1,  :-2, -1 ])
                                        / (2*DY)  +
                (w[ 1:-1, 1:-1,  0 ] - w[ 1:-1, 1:-1, -2 ])
                                        / (2*DZ)
        )

                                )
        ##set first equal to last
    src[ :, :, 0 ] = src[ :, :, -1 ]
    
    return src
# end set_bnds_periodic_src_div_only_3dZ
