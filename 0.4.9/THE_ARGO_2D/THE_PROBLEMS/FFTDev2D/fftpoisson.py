import numpy as np

def fftpoisson(src,dx):
    interior = src[1:-1,1:-1]
    [Nx,Ny] = interior.shape
    kx,ky = np.meshgrid(np.fft.fftfreq(Nx,d=dx),np.fft.fftfreq(Ny,d=dx))
    kx[0,0] = 1
    ky[0,0] = 1
    V = np.real_if_close(np.fft.ifft2(-np.fft.fft2(interior)/(kx**2+ky**2)))/(2*np.pi)**2
    return V    
    
