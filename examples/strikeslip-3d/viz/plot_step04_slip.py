#!/usr/bin/env python3

import h5py
import numpy as np
import matplotlib.pyplot as plt

# functions for plotting the fault mesh live here:
import plotFaultFunctions

def do_main(slipfname,gpsfname):

    with h5py.File(slipfname, 'r') as f:  # reads the hdf5 into a dictionary f
        vertices = np.array(f['geometry']['vertices'])
        cells = np.array(f['viz']['topology']['cells']) 
        slip = np.array(f['vertex_fields']['slip'])
    strike_slip = slip[:,:,1]
    reverse_slip = slip[:,:,2]
    slip_mag = np.sqrt(strike_slip**2 + reverse_slip**2)

    # Station Observations
    with h5py.File(gpsfname, "r") as h5:
        gps_coords = np.array(h5['geometry/vertices']) # Nstations x Ncomponents
        gps_disp = np.array(h5['vertex_fields/displacement']).squeeze() # Nstations x Ncomponents

    print(np.shape(gps_coords))
    print(np.shape(gps_disp))
    # plot the slip model in 2D, looking from above, 
    # using triangles to help improve the color interpolation

    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,5))
    plotFaultFunctions.plot_fault_surface(vertices,cells,strike_slip,ax1,'Strike Slip')
    plotFaultFunctions.plot_fault_surface(vertices,cells,reverse_slip,ax2,'Reverse slip')
    plt.show()

    # plot the predicted displacements
    fig,ax=plt.subplots()
    plotFaultFunctions.plot_fault_surface(vertices,cells,slip_mag,ax,'Slip Magnitude')
    ptsize=40
    # scatter plot for z values
    zscatter = ax.scatter(gps_coords[:,0], gps_coords[:,1], ptsize, gps_disp[:,2], cmap = 'RdBu_r',edgecolors='k')
    cbar = plt.colorbar(zscatter,ax=ax,extend='both')
    cbar.ax.set_title('Z (m)')
    # quiver for x and y values
    ax.quiver(gps_coords[:,0], gps_coords[:,1], gps_disp[:,0], gps_disp[:,1])
    
    plt.show()

if __name__ == "__main__":
    # enter the number of vertices in the mesh, in the x and y directions 
    slipfname="output/step04_varslip-fault.h5"
    gpsfname="output/step04_varslip-gps_stations.h5"

    do_main(slipfname,gpsfname)
    