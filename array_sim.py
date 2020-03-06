import numpy as np
import math, sys
import matplotlib.pyplot as plt


def antenna1(N, M, dx, dy, asrec = False, plot=False, show=False):
    
    """Function to create a mesh array of coordinates X,Y,Z
        
        Parameters
        ----------
        N: number of antennas along the x-axis
        M: number of antennas along the y-axis
        dx: antenna spacing in the x-axis
        dy: antenna spacing along the y-axis
        
        Returns
        --------
        - a list of antenna coordinates of X,Y,Z and ID if asrec is True.
        - Also return a mesh grid of the X, Y coordinates
        """
    
    x = dx*np.arange(0, N, dtype = 'float64')
    y = dy*np.arange(0, M, dtype = 'float64')
    
    
    dt = [('Pnt_num', '<i4'), ('X', '<f8'), ('Y', '<f8'),  ('Z', '<f8')]
    mesh = np.meshgrid(x, y, sparse = False)
    
    if asrec:
        xx = mesh[0].ravel()
        yy = mesh[1].ravel()
        zz = np.zeros(len(xx))
        
        arr = list(zip(np.arange(len(xx)), xx, yy, zz))
        pos_ants = np.array(arr, dtype=dt)
    else:
        arr = list(zip(mesh[0].ravel(), mesh[1].ravel(), np.zeros(len(mesh[0].ravel()))))
        pos_ants = np.array(arr)

    if plot:
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.scatter(pos_ants[:,0], pos_ants[:,1], color='b', linewidth=1,
                   linestyle='-', alpha=1);
        plt.axis('tight'); plt.grid(False);
        plt.xlabel(r'$x$', fontsize=16); plt.ylabel(r'$Y$', fontsize=16);
        if(show): plt.title(r'$\mathrm{Array}\ \mathrm{Layout}$', fontsize=18); plt.show()
    
    return pos_ants





def antenna_sim2(N,M, dx, dy, plot=False, show=False, write=True):
    
    """function generating array starting at the center
        
        Parameters
        ----------
        N: number of antennas along the x-axis
        M: number of antennas along the y-axis
        dx: antenna spacing in the x-axis
        dy: antenna spacing along the y-axis
        """
    [X,Y] = np.mgrid[0:N, 0:M]
    
    X = (X - (N - 1) / 2).T
    Y = np.flipud((Y - (N - 1) / 2).T)
    X = X * dx
    Y = Y * dy
    #arranging antenna position into [N, 3] matrix
    X = X.flatten()[:, None]
    Y = Y.flatten()[:, None]
    Z = np.zeros_like(X)
    
    xyz_pos = np.hstack((X,Y,Z))
    
    if write:
        header = "XYZ"
        with open('HIRAX_Antenna_Positions.txt', 'w+') as datafile_id:
            np.savetxt(datafile_id, xyz_pos, fmt=['%d','%d', '%d'], header = header)
    if plot:
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.scatter(X, Y, color='b', linewidth=1,
        linestyle='-', alpha=1);
        plt.axis('tight'); plt.grid(False),
        plt.xlabel(r'$x$', fontsize=16); plt.ylabel(r'$Y$', fontsize=16)
        if(show): plt.title(r'$\mathrm{Array}\ \mathrm{Layout}$', fontsize=18); plt.show()

    return xyz_pos

