#!/usr/bin/env python3

import h5py
import numpy as np
import matplotlib.pyplot as plt

# functions for plotting the fault mesh live here:
import plotFaultFunctions

def do_main(filename):

    with h5py.File(filename, 'r') as f:  # reads the hdf5 into a dictionary f
        vertices = np.array(f['geometry']['vertices'])
        cells = np.array(f['viz']['topology']['cells']) 
        slip = np.array(f['vertex_fields']['slip'])
    
    strike_slip = slip[:,:,1]
    reverse_slip = slip[:,:,2]
    # plot the slip model in 2D, looking from above, 
    # using triangles to help improve the color interpolation

    # Create the plot
    fig, ax1, ax2 = plt.subplots(1, 2,figsize=(12,5))
    plotFaultFunctions.plot_fault_surface(vertices,cells,strike_slip,ax1,'Strike Slip')
    plotFaultFunctions.plot_fault_surface(vertices,cells,reverse_slip,ax2,'Reverse slip')
    plt.show()

if __name__ == "__main__":
    # enter the number of vertices in the mesh, in the x and y directions 
    slipfname="output/step04_varslip-fault.h5"
    do_main(slipfname)
    