## burgers_1d.py
## by Ryan Farber 28 December 2014
## Last modified: 22 January 2015
"""
The purpose of this program is to apply the_argo to propagate a hat
function by burgers equation (nonlinear convection and diffusion)
in one dimension.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, "../../")
import the_bnds_setter_1d as bs 
sys.path.insert(0,"../../..")
from the_file_name_getter   import get_file_name
from the_state_saver        import save_state
from the_fluid              import The_Fluid
import os; cur_path = os.getcwd(); os.chdir(cur_path + "/StateFiles")
import cPickle; my_fluid = cPickle.load(open("my_fluid.p", "rb"))
X = np.linspace(my_fluid.XMIN, my_fluid.XMAX, my_fluid.NX)
NU = my_fluid.NU



##plot initial u_analytic
time = 0.0
phi = (np.exp(-(X - 4*time)**2 / (4*NU*(time+1)))
                +  np.exp(-(X - 4*time - 2*np.pi)**2 / (4*NU*(time+1))))

dphi_dx = (-X / (2*NU*(time+1)) * np.exp(-(X - 4*time)**2
                                         / (4*NU*(time+1)))
        - (X - 4*time - 2*np.pi) / (2*NU*(time+1)) * np.exp(
                    -(X - 4*time - 2*np.pi)**2 / (4*NU*(time+1))))

u_analytic = -(2*NU/phi)*dphi_dx + 4.0


##Setup
ext = ".p"  # filename extension for pickling
u = u_analytic


##Plot initial condition of fluid
cycles = 1; file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)

save_state( [cycles, u], file_name )

plt_title = "Sawtooth at cycle " + str(cycles)
plt_name  = file_name[:-2] + ".png"
os.chdir(cur_path + "/Plots")
plt.plot(X, u, marker='o', lw=2, label='Computational')
plt.plot(X, u_analytic, label='Analytical')
plt.xlim(my_fluid.XMIN, my_fluid.XMAX)
plt.ylim(0, 10)
plt.title(plt_title)
plt.savefig(plt_name)

plt.clf()


##Solve!
for cycles in xrange(1, my_fluid.NT+1):
    u = my_fluid.nonlinear_convect_1d(u)
    u = my_fluid.diffuse_1d(u)
    u = bs.set_bnds_periodic_diffusion_1d(my_fluid.DT, my_fluid.DX,
                                          my_fluid.NU, u)
    u = bs.set_bnds_fixed_1dI(u, u[-1])
# end for


##plot u_analytic
time = my_fluid.NT * my_fluid.DT
phi = (np.exp(-(X - 4*time)**2 / (4*NU*(time+1)))
                +  np.exp(-(X - 4*time - 2*np.pi)**2 / (4*NU*(time+1))))

dphi_dx = (-(X - 4*time) / (2*NU*(time+1)) * np.exp(-(X - 4*time)**2
                                         / (4*NU*(time+1)))
        - (X - 4*time - 2*np.pi) / (2*NU*(time+1)) * np.exp(
                -(X - 4*time - 2*np.pi)**2 / (4*NU*(time+1))))

u_analytic = -(2*NU/phi)*dphi_dx + 4.0


##plot computational result
file_name = get_file_name(cycles, my_fluid.MAX_CYC_MAG,
                                      my_fluid.LABEL, ext)

os.chdir(cur_path + "/StateFiles")
save_state( [cycles, u], file_name )

plt_title = "Sawtooth at cycle " + str(cycles)
plt_name  = file_name[:-2] + ".png"
os.chdir(cur_path + "/Plots")
plt.plot(X, u, marker='o', lw=2, label='Computational')
plt.plot(X, u_analytic, label='Analytical')
plt.xlim(my_fluid.XMIN, my_fluid.XMAX)
plt.ylim(0, 10)
plt.title(plt_title)
plt.savefig(plt_name)


## end burgers_1d.py
