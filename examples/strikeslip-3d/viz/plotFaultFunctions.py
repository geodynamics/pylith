#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def plot_fault_surface(vertices,cells,slipvec,ax,title):
    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, 2]
    slipvec=slipvec.flatten()
    
    # check if this is a quad or tri mesh surface - if quad, convert to triangles
    if np.shape(cells)[1] == 4:
        # List to store triangles
        triangles = []
        # Go through the cells and create triangles
        for cell in cells:
            idx0, idx1, idx2, idx3 = cell
            # Append center point of each cell to the lists
            x=np.append(x,np.mean(x[cell]))
            y=np.append(y,np.mean(y[cell]))
            z=np.append(z,np.mean(z[cell]))
            slipvec=np.append(slipvec,np.mean(slipvec[cell]))
            idx4 = len(x)-1
            # Divide each cell into four triangles
            triangles.append([idx0, idx1, idx4])
            triangles.append([idx1, idx2, idx4])
            triangles.append([idx2, idx3, idx4])
            triangles.append([idx3, idx0, idx4])
    else:
        # for a triangular mesh
        triangles=cells
    
    # create the image
    surf = ax.tricontourf(x,y,triangles,slipvec,20,extend='both')
    # Add a color bar which maps values to colors
    cbar = plt.colorbar(surf, ax=ax)
    # Add labels
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    cbar.ax.set_title('Slip (m)')
    ax.set_title(title)
    # plot lines showing the triangles, if desired
    ax.triplot(x,y,triangles,'w-',linewidth=0.1)
    