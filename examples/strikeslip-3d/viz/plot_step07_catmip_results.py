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
FILENAME_greensfns = "output/step05_greensfns-gps_stations.h5"
FILENAME_observed = "output/step04_varslip-gps_stations.h5"
CATMIP_results = "output/step07a_catmip-theta90.bin"
show_plot=True # change to False if running without graphics

## Read the appropriate files

# Station Observations
with h5py.File(FILENAME_observed, "r") as h5:
    observed_coords = h5['geometry/vertices'][:] # Nstations x Ncomponents
    observed_displacement = h5['vertex_fields/displacement'][:,:,:].squeeze() # Nstations x Ncomponents

# Greens functions
with h5py.File(FILENAME_greensfns, 'r') as h5:  # reads the hdf5 into a dictionary f
    greens_fns = np.array(h5['vertex_fields']['displacement'])

# un-interleave the greens functions
# could add a rotation matrix here too
greens_fns_reordered = np.vstack((greens_fns[0::2,:,:],greens_fns[1::2,:,:]))

# Impulse locations and components
with h5py.File(FILENAME_impulses, "r") as h5:
    impulse_coords = h5['geometry/vertices'][:]
    impulse_slip = h5['vertex_fields/slip'][:,:,:].squeeze() # Nimpulses x Ntotalnodes x Ncomponents
    vertices = np.array(h5['geometry']['vertices'])
    cells = np.array(h5['viz']['topology']['cells'])
    nimpulses = impulse_slip.shape[0]

# un-interleave the impulses
# could add a rotation matrix here too
impulse_slip_reordered = np.vstack((impulse_slip[0::2,:,:],impulse_slip[1::2,:,:]))

# CATMIP results output
with open(CATMIP_results, "rb") as fin:
    coefs = np.frombuffer(fin.read(), dtype=np.float64)
inversion_coefs = coefs.reshape((nimpulses,-1))

# median
median_coefs = np.median(inversion_coefs, axis=1)

# 95 % confidence interval
conf_interval_coefs =  np.percentile(inversion_coefs, 97.5, axis=1) - np.percentile(inversion_coefs, 2.5, axis=1)

# ordering coefs to all nodes in case greens functions are calculated only on a subset of available nodes
node_median_coefs = np.tensordot(median_coefs, impulse_slip_reordered, 1)
node_conf_interval_coefs = np.tensordot(conf_interval_coefs, impulse_slip_reordered, 1)

# Plot median and confidence interval for each component
fig, ax = plt.subplots(2, 2, figsize=(12,9))
plotFaultFunctions.plot_fault_surface(vertices,cells,node_median_coefs[:,1],ax[0,0],'Slip component 1')
plotFaultFunctions.plot_fault_surface(vertices,cells,node_median_coefs[:,2],ax[0,1],'Slip component 2')
plotFaultFunctions.plot_fault_surface(vertices,cells,node_conf_interval_coefs[:,1],ax[1,0],'$95\%$ CI Slip component 1')
plotFaultFunctions.plot_fault_surface(vertices,cells,node_conf_interval_coefs[:,2],ax[1,1],'$95\%$ CI Slip component 2')

figureFilename = "output/step07-catmip-results.pdf"
fig.savefig(figureFilename)

if show_plot:
    plt.show()
else:
    plt.close()

# compute the predicted displacements
median_pred_disp = np.tensordot(median_coefs, greens_fns_reordered, 1)

# get slip magnitude for plotting
slip_mag = np.sqrt(node_median_coefs[:,1]**2 + node_median_coefs[:,2]**2)

# plot the gps data
ptsize1=120
ptsize2=40
# scatter plot - observed z values
fig,ax=plt.subplots()
plotFaultFunctions.plot_fault_surface(vertices,cells,slip_mag,ax,'Slip magnitude')

zscatter = plt.scatter(observed_coords[:,0], observed_coords[:,1], ptsize1, observed_displacement[:,2], cmap = 'RdBu_r',edgecolors='k')
# scatter plot - predicted z values
zscatter = plt.scatter(observed_coords[:,0], observed_coords[:,1], ptsize2, median_pred_disp[:,2], cmap = 'RdBu_r',edgecolors='k')
# colorbar for scatter plot
cbar = plt.colorbar(zscatter,extend='both')
cbar.ax.set_title('Z (m)')
# quiver plot - observed displacements
plt.quiver(observed_coords[:,0], observed_coords[:,1], observed_displacement[:,0], observed_displacement[:,1])
# quiver plot - predicted displacements
plt.quiver(observed_coords[:,0], observed_coords[:,1], median_pred_disp[:,0], median_pred_disp[:,1],color='r')

figureFilename = "output/step07-catmip-gps_vectors.pdf"
fig.savefig(figureFilename)

if show_plot:
    plt.show()
else:
    plt.close()

# End of file
