# Step 2: Interseismic Deformation Simulation

In this example we simulate the interseismic deformation associated with the oceanic crust subducting beneath the continental crust and into the mantle.
We prescribe steady aseismic slip of 8 cm/yr along the interfaces between the oceanic crust and mantle with the interface between the oceanic crust and continental crust locked as shown in {numref}`fig:example:subduction:2d:steps`.
We adjust the Dirichlet boundary conditions on the lateral edges and bottom of the domain by pinning only the portions of the boundaries in the mantle and continental crust (i.e., not part of the oceanic crust).
Parameter settings that augment those in `pylithapp.cfg` are contained in the file `step02.cfg`.
These settings include:

**pylithapp.timedependent.formulation.time_step** Adjust the total simulation time to 100 years.

**pylithapp.timedependent** Specifies the array of boundary conditions.

**pylithapp.timedependent.bc.*BOUNDARY*** Defines the settings for boundary *BOUNDARY*, including which degrees of freedom are being constrained (x or y), the label (defined in` mesh_tri3.exo`) corresponding to the nodeset in CUBIT, and a label to the boundary condition used in any error messages.

**pylithapp.timedependent.interfaces** Specify the steady aseismic slip as a constant slip rate on the fault surfaces.

**pylithapp.problem.formulation.output.domain** Gives the base filename for HDF5 output (for example, `step02.h5`).

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02.cfg
```

The simulation will produce pairs of HDF5/Xdmf files with separate files for each material and fault interface.
{numref}`fig:example:subduction:2d:step02`, which was created using ParaView, displays the magnitude of the displacement field with the deformation exaggerated by a factor of 1000.
Using the animation features within ParaView or Visit you can illustrate how the continental crust near the trench subsides during the interseismic deformation.

:::{figure-md} fig:example:subduction:2d:step02
<img src="figs/step02_soln.*" alt="Solution for Step 2 at 100 years. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 2 at 100 years. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000.
:::
