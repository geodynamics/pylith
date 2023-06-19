#!/usr/bin/env python3

import h5py
import numpy as np
import matplotlib.pyplot as plt

# functions for plotting the fault mesh live here:
import plotFaultFunctions

def do_main(gpsfname,slipfname,chosen_idx):
    # Create the plot
    fig=plt.figure()
    ax=fig.subplots()

    # read slip data
    with h5py.File(slipfname, 'r') as f:  # reads the hdf5 into a dictionary f
        vertices = np.array(f['geometry']['vertices'])
        cells = np.array(f['viz']['topology']['cells']) 
        slip = np.array(f['vertex_fields']['slip'])
    strike_slip = slip[chosen_idx,:,1]
    reverse_slip = slip[chosen_idx,:,2]
    slip_mag=np.sqrt(strike_slip**2 + reverse_slip**2)
    # plot the slip model in 2D, looking from above, 
    # using triangles to help improve the color interpolation
    plotFaultFunctions.plot_fault_surface(vertices,cells,slip_mag,ax,'Greens Function %d'%chosen_idx)

    # read gps data
    with h5py.File(gpsfname, 'r') as f:  # reads the hdf5 into a dictionary f
        station_vertices = np.array(f['geometry']['vertices']);  # shape = (50, 3)
        gps_disp = np.array(f['vertex_fields']['displacement'])  # shape = (200, 50, 3)
    # View a given patch's greens function
    station_geo_x = station_vertices[:, 0]
    station_geo_y = station_vertices[:, 1]
    u_disp = gps_disp[chosen_idx, :, 0]
    v_disp = gps_disp[chosen_idx, :, 1]
    w_disp = gps_disp[chosen_idx, :, 2]
    
    ptsize=40
    # scatter plot for z values
    zscatter = plt.scatter(station_geo_x, station_geo_y, ptsize, w_disp, cmap = 'PuOr_r',edgecolors='k')
    cbar = plt.colorbar(zscatter,extend='both')
    cbar.ax.set_title('Z (m)')
    # quiver for x and y values
    plt.quiver(station_geo_x, station_geo_y, u_disp, v_disp)
    plt.show()

if __name__ == "__main__":
    
    # choose which greens fn to plot
    chosen_idx = 130;

    gpsfname="output/step05_greensfns-gps_stations.h5"
    slipfname="output/step05_greensfns-fault.h5"
    do_main(gpsfname,slipfname,chosen_idx)
    