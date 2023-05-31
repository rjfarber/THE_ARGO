import numpy as np

def jacobipoisson(src,dx):
    Nx,Ny = src.shape
    V = np.zeros((Nx,Ny))
    for i in xrange(0,int(1e4)+1):
        V[1:-1,1:-1] = 0.1*V[1:-1,1:-1] + 0.9*(V[:-2,1:-1]+V[2:,1:-1]+V[1:-1,:-2]+V[1:-1,2:] - src[1:-1,1:-1]*dx**2) / 4.0
        V[:,0] = V[:,-2]; V[:,-1] = V[:,1]
        V[0,:] = V[-2,:]; V[-1,:] = V[1,:]
    # end
    return V[1:-1,1:-1]
