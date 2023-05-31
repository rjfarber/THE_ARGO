## jacobi_vs_fft_2d.py
## by Ryan Farber 19 March 2015
## Last modified: 04 April 2015
"""
The purpose of this program is to compare the_argo's different pressure
solvers for a top hat src function in two dimensions.
"""
import numpy as np
import sys; sys.path.insert(0, "../../..")    
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; os.chdir("./StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))

##Setup
ext = ".p"  # filename extension for pickling
p   = np.zeros( (my_fluid.NX,my_fluid.NY) ) # pressure
src = np.zeros( (my_fluid.NX,my_fluid.NY) ) # source term of poisson eqn

src[0.25*my_fluid.NX : 0.75*(my_fluid.NY+1),
    0.25*my_fluid.NX : 0.75*(my_fluid.NY+1)] = 1


##Save the state of the initial condition of the fluid
cycles = 0; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                               my_fluid.LABEL, ext)
save_state([cycles, "NA", "NA","NA", p,src, "NA","NA", "NA","NA"],
           file_name )

##Jacobi
pct_err = 1.0; cycles += 1
for i in xrange(0,int(1e4)+1):
    p[1:-1,1:-1] = 0.1*p[1:-1,1:-1] + 0.9*(p[ :-2,1:-1]+p[2:  ,1:-1]
                                          +p[1:-1, :-2]+p[1:-1,2:  ]
            -src[1:-1,1:-1]*my_fluid.DX*my_fluid.DY) / 4.0
    
    #p = my_fluid.relax_pressure_poisson_2d(p, src)

    p[ 0, : ] = p[ -2,  : ]; p[ -1,  : ] = p[ 1, : ]
    p[ :, 0 ] = p[  :, -2 ]; p[  :, -1 ] = p[ :, 1 ]
# end for
pJac = p - np.mean(p)
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                   my_fluid.LABEL, ext)
save_state([cycles, "NA", "NA","NA", pJac,src, "NA","NA", "NA","NA"],
           "jacobi"+file_name )

##FFT
p = np.zeros( (my_fluid.NX,my_fluid.NY) )
p = my_fluid.transform_pressure_poisson_2d(p, src)
p[ 0, : ] = p[ -2,  : ]; p[ -1,  : ] = p[ 1, : ]
p[ :, 0 ] = p[  :, -2 ]; p[  :, -1 ] = p[ :, 1 ]
pFFT2 = p - np.mean(p)
save_state([cycles, "NA", "NA","NA", pFFT2,src, "NA","NA", "NA","NA"],
           "FFT2"+file_name )
print('Relative error:')
print("pFFT2 from pJac: ", np.max(np.absolute(pFFT2-pJac)/np.mean(np.mean(np.absolute(pJac)))))

## end jacobi_vs_fft_2d.py
