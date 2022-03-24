# Step 7: Inversion of Slow-Slip Event using 3-D Green's Functions

This example is a three-dimensional analog of {ref}`sec:example:greensfns2d` and is a more realistic example of how PyLith can be used to perform geodetic inversions.
We divide generating Green's functions for slip impulses on the central rupture patch of the subduction interface two sub-problems:

**Step 7a**  Left-lateral slip component.

**Step 7b**  Reverse slip component.

Although PyLith can generate the two components in one simulation, we often prefer to speed up the process by running simulations for each of the components at the same time using multiple processes on a cluster.

To generate the Green's functions we change the problem from the default *TimeDependent* to *GreensFns*.
We do this on the command line (as illustrated below).
PyLith automatically reads the `greensfns.cfg` parameter file.
This file contains settings that are common to both sub-problems.
 Note that the settings in the `greensfns.cfg` only apply to parameters associated with the *GreensFns* and its sub-components.
 For the Green's function problem, we must specify the fault interface and the id for the fault.
 We specify the amplitude of the impulses via a *UniformDB* spatial database, because we want impulses over the entire fault patch.
 We also request the amplitude of the impulses to be included in the fault info file.

```{code-block} cfg
---
caption: Excerpt from `greensfns.cfg`
---
# Define the interfaces (slab) and provide a fault_id.
[greensfns]
interfaces = [slab]
fault_id = 100

# Switch fault to FaultCohesiveImpulses for generation of Green's functions.
[greensfns.interfaces]
slab = pylith.faults.FaultCohesiveImpulses

[greensfns.interfaces.slab]
# Nodesets corresponding to the fault and its buried edges.
label = fault_slabtop_patch
edge = fault_slabtop_patch_edge

# We must define the quadrature information for fault cells.
# The fault cells are 2D (surface).
quadrature.cell = pylith.feassemble.FIATSimplex
quadrature.cell.dimension = 2

# Spatial database for slip impulse amplitude.
db_impulse_amplitude = spatialdata.spatialdb.UniformDB
db_impulse_amplitude.label = Amplitude of fault slip impulses
db_impulse_amplitude.values = [slip]
db_impulse_amplitude.data = [1.0]

# Add impulse amplitude to fault info output.
output.vertex_info_fields = [normal_dir, strike_dir, dip_dir, impulse_amplitude]
output.writer = pylith.meshio.DataWriterHDF5
```

We do not make use of the state variable output for the impulse responses, so we turn off the data fields for all of the materials to eliminate these large data files.

```{code-block} cfg
---
caption: Excerpt from `greensfns.cfg`
---
# Turn off output of state variables for materials.
[greensfns.materials.slab.output]
cell_data_fields = []

[greensfns.materials.wedge.output]
cell_data_fields = []

[greensfns.materials.crust.output]
cell_data_fields = []

[greensfns.materials.mantle.output]
cell_data_fields = []
```

The `step07a.cfg` and `step07b.cfg` files are identical, except for the impulse type specification and file names.

```{code-block} cfg
---
caption: Excerpt from `step07a.cfg`
---
[pylithapp.problem.interfaces.slab]
# If we wanted to generate impulses for both the left-lateral and
# reverse components in the same simulation, we would use:
# impulse_dof = [0,1]
#
# Impulses for left-lateral slip.
impulse_dof = [0]
```

In the output settings, we turn off writing the solution field for the domain:

```{code-block} cfg
---
caption: Excerpt from `step07a.cfg`
---
[pylithapp.problem.formulation.output.domain]
writer.filename = output/step07a-domain.h5
# Turn off data fields.
vertex_data_fields = []
```

```{code-block} console
---
caption: Run Step 7 simulations
---
$ pylith --problem=pylith.problems.GreensFns step07a.cfg mat_elastic.cfg solver_fieldsplit.cfg
$ pylith --problem=pylith.problems.GreensFns step07b.cfg mat_elastic.cfg solver_fieldsplit.cfg
```
Each simulation will produce four pairs of HDF5/Xdmf files.
For Step 7a these will be:

**step07a-groundsurf.h5[.xmf]**  Solution field over the ground surface for each slip impulse.

**step07a-cgps_sites.h5[.xmf]**  Solution field at continuous GPS sites for each slip impulse.

**step07a-fault-slab_info.h5[.xmf]**  Fault orientation and impulse information.

**step07a-fault-slab.h5[.xmf]**  Fault slip for each slip impulse.

:::{tip}
To save time, run the two sub-problems simultaneously in separate shells (terminal windows or tabs).
For a problem this size, this should work fine on a laptop.
For larger problems, we would run the simulations via separate jobs on a cluster with each job running on multiple processes.
:::

Before we can run the inversion, we post-process the output from Step 6 to create synthetic data.
We use the same generalized inverse approach described in {ref}`sec:example:greensfns2d:inversion`.
The Python script `make_synthetic_gpsdisp.py` reads the parameters in `make_synthetic_gpsdisp.cfg` and generates synthetic data from the selected time step with a specified amount of noise.

```{code-block} console
---
caption: Generate synthetic GPS data
---
./make_synthetic_gpsdisp.py
```

This will create the following files:

**cgps_synthetic_displacement.txt** read by the inversion script.

**cgps_synthetic_displacement.vtk** for visualization.

We perform a simple inversion using the `slip_invert.py` script, with parameters defined in `slip_invert.cfg`.
This script performs a set of linear inversions, in a manner similar to the inversion in {ref}`sec:example:greensfns2d:inversion`.

```{code-block} console
$./slip_invert.py
```

This will create a number of files in the output directory.

**step07-inversion-slip.h5**  This HDF5/Xdmf pair of files may be used to visualize the predicted slip distributions for different values of the penalty weight.

**step07-inversion-displacement.h5**  This HDF5/Xdmf pair of files may be used to visualize the predicted cGPS displacements for each solution.

**step07-inversion-summary.txt**  This file provides a summary of the inversion results for each value of the penalty weight.


One approach to finding the optimal penalty weight is to find the corner of the 'L-curve' for the log of the weighted data residual versus the log of the penalty residual.
This is viewed as the point of diminishing returns for reducing the penalty weight.
Further reductions provide little improvement to the weighted data residual, while providing a solution with less regularization.
{numref}`fig:example:subduction:3d:step07:curve` shows that this procedure suggests an optimal penalty weight of 0.1 for our inversion.

:::{figure-md} fig:example:subduction:3d:step07:curve
<img src="figs/subduction3d_step07_inverse_curve.*" alt="Plot of the 'L-curve' for inversion in Step 7. The 'corner' of the L-curve would be about the third or fourth point from the right of the plot, representing a penalty weight of 0.5 or 1.0 in our example." width="100%"/>

Plot of the 'L-curve' for inversion in Step 7. The 'corner' of the L-curve would be about the third or fourth point from the right of the plot, representing a penalty weight of 0.5 or 1.0 in our example.
:::

{numref}`fig:example:subduction:3d:step07:soln` shows the predicted slip, the observed and predicted displacement vectors, and the slip applied from example step06 for a penalty weight of 1.0.
The data fit is very good, and the predicted slip distribution is very close to the applied slip, although the magnitude is slightly underestimated.

:::{figure-md} fig:example:subduction:3d:step07:soln
 <img src="figs/subduction3d_step07_inverse_soln.*" alt="ParaView image of the inversion solution for a penalty weight of 1.0. 'Data' is shown with blue arrows and predicted displacements are shown with magenta arrows. Color contours represent the predicted slip distribution and orange line contours show the applied slip from the forward problem." width="100%"/>

ParaView image of the inversion solution for a penalty weight of 1.0. 'Data' is shown with blue arrows and predicted displacements are shown with magenta arrows. Color contours represent the predicted slip distribution and orange line contours show the applied slip from the forward problem.
:::

## Exercises

* Investigate the effects of data noise.
  * How do the noisy data vectors compare to the raw data vectors from example step06?
  * Create a new simulated dataset with more noise and see how well the solution matches the applied slip.
* Different initial slip distribution.
  * Move the slip distribution to a different location, vary the amplitude, etc. This will involve running another instance of example step06 to create a new dataset. How is the solution affected?
  * Move the slip onto the splay fault. This will involve creating a new forward model as well as generating Green's functions for the splay fault.
* What happens if your material properties are incorrect?
  * Try creating your forward model with heterogeneous properties and your Green's functions with homogeneous properties (or vice-versa). What happens to your solution?
* Try inverting for slip at various time steps.
* Try a different inversion method.
  * If you analyze the predicted slip distribution you will find some negative slip, which is unrealistic. To overcome this problem you could try NNLS (non-negative least squares). If you have the Python scipy package installed on your computer, you could replace the generalized inverse solution with the NNLS package included in `scipy.optimize.nnls`.
