import numpy as np

def fftpoisson2(src,dx):
    interior = src[1:-1,1:-1]
    [Nx,Ny] = interior.shape
    kx,ky = np.meshgrid(np.fft.fftfreq(Nx,d=dx),np.fft.fftfreq(Ny,d=dx))
    factor = 1.0/(4 - 2*np.cos(2*np.pi*kx*dx) - 2*np.cos(2*np.pi*ky*dx))
    factor[0,0] = 0
    V = np.real_if_close(np.fft.ifft2(-dx**2 * np.fft.fft2(interior)*factor))
    return V    
    
