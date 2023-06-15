"""
This script generates a plot showing the median and 95% confidence interval of the posterior samples from catmip
Authors: Rishav Mallick & Eric Lindsey
"""

# Import numpy, h5py, and matplotlib modules (included in PyLith binary installation)
import numpy as np
import h5py
import matplotlib.pyplot as plt

# functions for plotting the fault mesh live here:
import plotFaultFunctions

FILENAME_impulses = "output/step05_greensfns-fault.h5"
FILENAME_observed = "output/step04_varslip-gps_stations.h5"
CATMIP_results = "output/step07a_catmip-theta45.bin"
show_plot=True # change to False if running without graphics

## Read the appropriate files

# Station Observations (currently unused) 
h5 = h5py.File(FILENAME_observed, "r")
observed_coords = h5['geometry/vertices'][:] # Nstations x Ncomponents
observed_displacement = h5['vertex_fields/displacement'][:,:,:].squeeze() # Nstations x Ncomponents
h5.close()

# Greens functions (impulse sources)
h5 = h5py.File(FILENAME_impulses, "r")
impulse_coords = h5['geometry/vertices'][:]
impulse_slip = h5['vertex_fields/slip'][:,:,:].squeeze() # Nimpulses x Ntotalnodes x Ncomponents
vertices = np.array(h5['geometry']['vertices'])
cells = np.array(h5['viz']['topology']['cells'])
h5.close()
nimpulses = impulse_slip.shape[0]

# un-interleave the greens functions
# could add greens fns rotation matrix here too
impulse_slip_reordered = np.vstack((impulse_slip[0::2,:,:],impulse_slip[1::2,:,:]))
print(np.shape(impulse_slip_reordered))

# CATMIP results output
with open(CATMIP_results, "rb") as fin:
    coefs = np.frombuffer(fin.read(), dtype=np.float64)
inversion_coefs = coefs.reshape((nimpulses,-1))

# median
median_coefs = np.median(inversion_coefs, axis=1)
print(np.shape(median_coefs))

# 95 % confidence interval
conf_interval_coefs =  np.percentile(inversion_coefs, 97.5, axis=1) - np.percentile(inversion_coefs, 2.5, axis=1)

# ordering coefs to all nodes in case greens functions are calculated only on a subset of available nodes
node_median_coefs = np.tensordot(median_coefs, impulse_slip_reordered, 1)
node_conf_interval_coefs = np.tensordot(conf_interval_coefs, impulse_slip_reordered, 1)

# Create the plot 
fig, ax = plt.subplots(2, 2, figsize=(12,9))
plotFaultFunctions.plot_fault_surface(vertices,cells,node_median_coefs[:,1],ax[0,0],'Slip component 1')
plotFaultFunctions.plot_fault_surface(vertices,cells,node_median_coefs[:,2],ax[0,1],'Slip component 2')
plotFaultFunctions.plot_fault_surface(vertices,cells,node_conf_interval_coefs[:,1],ax[1,0],'$95\%$ CI Slip component 1')
plotFaultFunctions.plot_fault_surface(vertices,cells,node_conf_interval_coefs[:,2],ax[1,1],'$95\%$ CI Slip component 2')

figureFilename = "output/step07-catmip-results.pdf"
fig.savefig(figureFilename)

if show_plot:
    plt.show()

# End of file
