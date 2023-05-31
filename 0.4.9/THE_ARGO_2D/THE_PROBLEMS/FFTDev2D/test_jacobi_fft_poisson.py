import numpy as np
from jacobipoisson import *
from fftpoisson    import *
from fftpoisson2   import *

N = 20
dx = 2.0/(N-2)
x = np.arange(-dx, (N-2)*dx, dx) # does not include 1.0
y = x
xx,yy = np.meshgrid(x,y) # I've confirmed the order of output is same as matlab

src = np.zeros((N,N))
src[N/4.0:(3.0*N/4 + 1), N/4.0:(3.0*N/4 + 1)] = 1
#src = src - np.mean(np.mean(src))

# For sine wave
#kx = 4*np.pi
#ky = 6*np.pi
#src = np.sin(kx*xx) * np.sin(ky*yy)

Vjac  = jacobipoisson(src,dx)

Vfft  = fftpoisson( src,dx)
Vfft2 = fftpoisson2(src,dx) 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fig1 = plt.figure()
ax = fig1.gca(projection='3d')
Vjac = Vjac - np.mean(np.mean(Vjac))

ax.plot_surface(xx[1:,1:],yy[1:,1:],Vjac, rstride=1,cstride=1, cmap=cm.coolwarm,
                linewidth=0, antialiased=False)
ax.set_title('Jacobi')
plt.savefig("NumpyJacobi.png")
print(np.min(Vjac))

fig1 = plt.figure()
ax = fig1.gca(projection='3d')
Vfft = Vfft - np.mean(np.mean(Vfft))
ax.plot_surface(xx[1:,1:],yy[1:,1:],Vfft, rstride=1,cstride=1, cmap=cm.coolwarm,
                linewidth=0, antialiased=False)
ax.set_title('FFT')
plt.savefig("NumpyFFT.png")
print(np.min(Vfft))

fig1 = plt.figure()
ax = fig1.gca(projection='3d')
Vfft2 = Vfft2 - np.mean(np.mean(Vfft2))
ax.plot_surface(xx[1:,1:],yy[1:,1:],Vfft2, rstride=1,cstride=1, cmap=cm.coolwarm,
                linewidth=0, antialiased=False)
ax.set_title('FFT2')
plt.savefig("NumpyFFT2.png")
print(np.min(Vfft2))

print('Relative error:')
print(" Vfft from Vjac: ", np.max(np.absolute(Vfft-Vjac))/np.mean(np.mean(np.absolute(Vjac))))
print("Vfft2 from Vjac: ", np.max(np.absolute(Vfft2-Vjac))/np.mean(np.mean(np.absolute(Vjac))))


