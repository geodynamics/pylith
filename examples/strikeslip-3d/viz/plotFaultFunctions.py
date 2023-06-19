#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def plot_fault_surface(vertices,cells,slipvec,ax,title):
    # INPUTS are 
    # vertices - Nx3 matrix of elements on the fault (x,y,z)
    # cells    - index array describing the connection of vertices to form fault polygons
    #            (can be Mx3 for triangles or Mx4 for quadrilaterals)
    # slipvec  - Nx1 vector of slip to use for coloring the fault
    # ax       - matplotlib fault axes object
    # title    - title to use for the axes

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
    vmax = np.ceil(np.max(np.abs(slipvec)))
    vmin = -vmax
    if vmin == vmax:
        vmin=-1
        vmax=1
    levels = np.linspace(vmin, vmax, 201)
    print(levels)
    surf = ax.tricontourf(x,y,triangles,slipvec,levels=levels,extend='both',cmap='RdBu_r')

    # Add a color bar which maps values to colors
    cbar = plt.colorbar(surf, ax=ax)
    # Add labels
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    cbar.ax.set_title('Slip (m)')
    ax.set_title(title)
    # plot lines showing the triangles, if desired
    ax.triplot(x,y,triangles,'w-',linewidth=0.1)

def computeUnitVectors(vertices,cells):
    # INPUTS are 
    # vertices - Nx3 matrix of elements on the fault (x,y,z)
    # cells - index array describing the connection of vertices to form fault polygons
    #         (can be Mx3 for triangles or Mx4 for quadrilaterals)

    # number of vertices and polygons
    nVertices=vertices.shape[0]
    nPatch=cells.shape[0]

    # vertex positions
    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, 2]

    # Go through the cells and create triangles
    i = 0
    nv = np.zeros((nPatch,3))
    sv = np.zeros((nPatch,3))
    dv = np.zeros((nPatch,3))
    xc = np.zeros((nPatch,3))

    for cell in cells:
        idx0, idx1, idx2 = cell[0:3]
        # vertices
        A = np.array((x[idx0],y[idx0],z[idx0]))
        B = np.array((x[idx1],y[idx1],z[idx1]))
        C = np.array((x[idx2],y[idx2],z[idx2]))

        # center points of each polygon
        xc[i,:]=np.array((np.mean(x[cell]),np.mean(y[cell]),np.mean(z[cell])))

        # normal vector (B-A) x (C-A)
        nv[i,:]=np.cross(B-A,C-A)
        area=np.sqrt(nv[i,0]**2+nv[i,1]**2+nv[i,2]**2)
        nv[i,:]=nv[i,:]/area

        # strike-direction vector
        sv[i,:]=np.array((-np.sin(np.arctan2(nv[i,1],nv[i,0])),
                        np.cos(np.arctan2(nv[i,1],nv[i,0])),
                        0))
        
        # dip-direction vector
        dv[i,:]=np.array((nv[i,1]*sv[i,2]-nv[i,2]*sv[i,1],
                        nv[i,2]*sv[i,0]-nv[i,0]*sv[i,2],
                        nv[i,0]*sv[i,1]-nv[i,1]*sv[i,0]))

        i+=1
    return xc,nv,sv,dv 