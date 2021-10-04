# Step 1: Coseismic Slip Simulation

The first example problem is earthquake rupture involving coseismic slip along the interface between the subducting slab and the continental crust and uppermost portion of the mantle below the continental crust.
The spatial variation of slip comes from a cross-section of Gavin Hayes's finite-source model <https://earthquake.usgs.gov/earthquakes/eventpage/usc0001xgp/finite-fault>.
On the lateral and bottom boundaries of the domain, we fix the degrees of freedom perpendicular to the boundary as shown in {numref}`fig:example:subduction:2d:steps`.
Parameter settings that augment those in `pylithapp.cfg` are contained in the file `step01.cfg`.
These settings are:

**pylithapp.timedependent.formulation.time_step** Adjust the total simulation time to 0 years (static simulation).

**pylithapp.timedependent** Specifies the array of boundary conditions.

**pylithapp.timedependent.bc.*BOUNDARY*** Defines the settings for boundary *BOUNDARY*, including which degrees of freedom are being constrained (x or y), the label (defined in `mesh_tri3.exo`) corresponding to the nodeset in CUBIT, and a label to the boundary condition used in any error messages.

**pylithapp.timedependent.interfaces.fault** Specify the coseismic slip along the interface between the oceanic crust and continental crust with a small amount of slip penetrating into the upper mantle.

**pylithapp.problem.formulation.output.domain** Gives the base filenames for HDF5 output (for example, `step01.h5`).

```{code-block} console
---
caption: Run Step 1 simulation
---
$ pylith step01.cfg
```

The problem will produce twelve pairs of HDF5/Xdmf files.
The HDF5 files contain the data and the Xdmf files contain the metadata required by ParaView and Visit (and possibly other visualization tools that use Xdmf files) to access the mesh and data sets in the HDF5 files.
The files include the solution over the domain and ground surface (two pairs of files), physical properties, stress, and strain within each material (eight pairs of files), and fault parameters, slip, and traction (two pairs of files).

{numref}`fig:example:subduction:2d:step01`, which was created using ParaView, displays the magnitude of the displacement field with the deformation exaggerated by a factor of 1000.

:::{figure-md} fig:example:subduction:2d:step01
<img src="figs/step01_soln.*" alt="Solution for Step1. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000." width="100%">

Solution for Step 1. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000.
:::
