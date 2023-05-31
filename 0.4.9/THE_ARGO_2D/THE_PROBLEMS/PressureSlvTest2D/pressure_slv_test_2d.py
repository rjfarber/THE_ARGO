## pressure_slv_test_2d.py
## by Ryan Farber 05 March 2015
## Last modified: 18 March 2015
"""
The purpose of this program is to apply the_argo to propagate
the orszag-tang vortex implicitly, a 2D ideal MHD problem,
in incompressible form.

I am attempting to imitate as closely as possible:
https://perswww.kuleuven.be/~u0016541/MHD_sheets_pdf/AdvCFDles04.pdf
"""
import glob
import numpy as np
import sys; sys.path.insert(0, "../../..")    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Setup
ext = ".p"  # filename extension for pickling
X = np.linspace(my_fluid.XMIN,my_fluid.XMAX,my_fluid.NX)
Y = np.linspace(my_fluid.YMIN,my_fluid.YMAX,my_fluid.NY)
Y,X = np.meshgrid(X,Y)

XX,YY = np.mgrid[1:my_fluid.NX-1, 1:my_fluid.NY-1] #imitates 2d for loop

if my_fluid.cycle_start == 1:
    u   = np.ones( (my_fluid.NX,my_fluid.NY)) # x-component of velocity
    v   = np.zeros((my_fluid.NX,my_fluid.NY)) # y-component of velocity
    cn  = np.zeros((my_fluid.NX,my_fluid.NY)) # concentration (tracers)
    p   = np.zeros((my_fluid.NX,my_fluid.NY)) # pressure
    src = np.zeros((my_fluid.NX,my_fluid.NY)) # src term for poisson eqn
    Bx  = np.zeros((my_fluid.NX,my_fluid.NY)) # x-comp of magnetic field
    By  = np.zeros((my_fluid.NX,my_fluid.NY)) # y-comp of magnetic field

    u[0.5/my_fluid.DX : 1.0/my_fluid.DX+1,
      0.5/my_fluid.DY : 1.0/my_fluid.DY+1 ] = 2.0

    ##Update ghost zones
    u[ 0, : ] = u[ -2,  : ]; u[ -1,  : ] = u[ 1, : ]
    u[ :, 0 ] = u[  :, -2 ]; u[  :, -1 ] = u[ :, 1 ]

    v[ 0, : ] = v[ -2,  : ]; v[ -1,  : ] = v[ 1, : ]
    v[ :, 0 ] = v[  :, -2 ]; v[  :, -1 ] = v[ :, 1 ]

    ##Save the state of the initial condition of the fluid
    cycles = 0; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                          my_fluid.LABEL, ext)
    save_state([cycles, u,v, cn, p,src, "NA","NA", Bx,By],
               file_name )
else:
    data_file = glob.glob("*" + str(my_fluid.cycle_start) + "*")
    if data_file == []:
        print("Error! my_fluid.cycle_start data file not found.")
        sys.exit()
    # end if
    data_file = data_file[0]
    the_data = cPickle.load(open(data_file))
    u  = the_data[1]; v   = the_data[2]; cn = the_data[3]
    p  = the_data[4]; src = the_data[5]; Bx = the_data[7]
    By = the_data[8]
    the_data = 0 # to save memory
# end if


##Solve!
for cycles in xrange(my_fluid.cycle_start, my_fluid.NT+1):
    divV  = (u[ 2:  , 1:-1 ] - u[  :-2, 1:-1 ]) \
                                      / (2*my_fluid.DX)
    divV += (v[ 1:-1, 2:   ] - v[ 1:-1,  :-2 ]) \
                                      / (2*my_fluid.DY)
    divV = sum(sum(np.absolute(divV)))

    div_iters = 0
    while divV > 1e-3:
        div_iters += 1
        src[1:-1,1:-1] = my_fluid.calc_source_2d(u,v)

        if div_iters == 1:
            print
            print("S: ", my_fluid.S)
            print("NI: ", my_fluid.NI)
            absMaxSrc = sum(sum(np.absolute(src)))
            print("before: ", absMaxSrc)
        
        for dummy_var in xrange(my_fluid.NI):        
            p = my_fluid.relax_pressure_poisson_2d(p, src)

            p[ 0, : ] = p[ -2,  : ]; p[ -1,  : ] = p[ 1, : ]
            p[ :, 0 ] = p[  :, -2 ]; p[  :, -1 ] = p[ :, 1 ]
        # end for

        u[1:-1,1:-1] -= my_fluid.apply_pressure_2dX(p, 1.0/my_fluid.RHO)
        v[1:-1,1:-1] -= my_fluid.apply_pressure_2dY(p, 1.0/my_fluid.RHO)

        ##Update ghost zones
        u[ 0, : ] = u[ -2,  : ]; u[ -1,  : ] = u[ 1, : ]
        u[ :, 0 ] = u[  :, -2 ]; u[  :, -1 ] = u[ :, 1 ]

        v[ 0, : ] = v[ -2,  : ]; v[ -1,  : ] = v[ 1, : ]
        v[ :, 0 ] = v[  :, -2 ]; v[  :, -1 ] = v[ :, 1 ]

        src[1:-1,1:-1] = my_fluid.calc_source_2d(u,v)
        absMaxSrc = sum(sum(np.absolute(src)))
        #print("after: ", absMaxSrc)
        #print
        
        divV  = (u[ 2:  , 1:-1 ] - u[  :-2, 1:-1 ]) \
                                      / (2*my_fluid.DX)
        divV += (v[ 1:-1, 2:   ] - v[ 1:-1,  :-2 ]) \
                                      / (2*my_fluid.DY)
        divV = sum(sum(np.absolute(divV)))
    # end while
    print("num div iters: ", div_iters)
    print("after: ", absMaxSrc)
    print

    if cycles == 4:
        sys.exit()


    if (cycles in my_fluid.SAVE_SP)or(cycles % my_fluid.SAVE_FREQ == 0):
        file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                           my_fluid.LABEL, ext)
        save_state([cycles, u,v, cn, p,src, "NA","NA", Bx,By],
                   file_name)
    # end if
# end for
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                   my_fluid.LABEL, ext)
save_state([cycles, u,v, cn, p,src, "NA","NA", Bx,By],
           file_name )

## end pressure_slv_test_2d.py
