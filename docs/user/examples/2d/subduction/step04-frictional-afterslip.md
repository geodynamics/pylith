# Step 4: Frictional Afterslip Simulation

This simulation demonstrates how to combine the change in tractions associated with coseismic slip with a background stress field to compute afterslip controlled by static friction.
The Python script `afterslip_tractions.py` will create a spatial database file with initial tractions based on the change in tractions from Step 1 and a background stress field.
The background stress field is simply normal tractions consistent with the overburden (lithostatic load) for a uniform half-space and shear tractions consistent with a coefficient of friction of 0.6.
The `afterslip_tractions.spatialdb` file is provided, so you do not need to run the Python script `afterslip_tractions.py`; however, you can do so by typing

```{code-block} console
---
caption: Generate `afterslip_tractions.spatialdb`
---
$ python afterslip_tractions.py
```

We provide 2.0 MPa of strength excess associated with the background stress field by using a cohesion of 2.0 MPa in the static friction model.
Slip will occur in regions where the coseismic slip increased the shear tractions by more than 2.0 MPa.
On the lateral and bottom boundaries of the domain, we fix the degrees of freedom perpendicular to the boundary as shown in {numref}`fig:example:subduction:2d:steps`.
Parameter settings that augment those in `pylithapp.cfg` are contained in the file `step04.cfg`. These settings are:

**pylithapp.timedependent.formulation.time_step** Adjust the total simulation time to 0 years (static simulation).

**pylithapp.timedependent** Selects the nonlinear solver and specifies the array of boundary conditions.

**pylithapp.timedependent.bc.*BOUNDARY*** Defines the settings boundary *BOUNDARY*, including which degrees of freedom are being constrained (x or y), the label (defined in `mesh_tri3.exo`) corresponding to the nodeset in CUBIT, and a label to the boundary condition used in any error messages.

**pylithapp.timedependent.interfaces** Specify a fault with a fault constitutive model (static friction) and initial fault tractions.

**pylithapp.problem.formulation.output.domain** Gives the base filename for HDF5 output (for example, `step04.h5`).

```{code-block} console
---
caption: Run Step 4 simulation
---
$ pylith step04.cfg
```

The problem will produce twelve pairs of HDF5/Xdmf files.
The HDF5 files contain the data and the Xdmf files contain the metadata required by ParaView and Visit (and possibly other visualization tools that use Xdmf files) to access the mesh and data sets in the HDF5 files.
The files include the solution over the domain and ground surface (two pairs of files), physical properties, stress, and strain within each material (eight pairs of files), and fault parameters, slip, and traction (two pairs of files).

{numref}`fig:example:subduction:2d:step04`, which was created using ParaView, displays the magnitude of the displacement field with the original configuration. Slip occurs down-dip from the coseismic slip as well as in three areas with sharp gradients in slip, including the trench.
The location of the afterslip can be shifted by changing the spatial variation of the coseismic slip and background stress field.

:::{figure-md} fig:example:subduction:2d:step04
<img src="figs/step04_soln.*" alt="Solution for Step 4. The colors indicate the magnitude of the displacement" width="100%"/>

Solution for Step 4. The colors indicate the magnitude of the displacement.
:::
