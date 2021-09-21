# Step 3: Pseudo-Earthquake Cycle Model

This simulation combines 300 years of interseismic deformation from Step 2 with the coseismic deformation from Step 1 applied at 150 years to create a simple model of the earthquake cycle.
Parameter settings that augment those in `pylithapp.cfg` are contained in the file `step03.cfg`.
These settings include:

**pylithapp.timedependent.formulation.time_step** Adjust the total simulation time to 300 years.

**pylithapp.timedependent** Specifies the array of boundary conditions.

**pylithapp.timedependent.bc.*BOUNDARY*** The Dirichlet boundary conditions match those in Step 2.

**pylithapp.timedependent.interfaces** On the interface between the subducting oceanic crust and the mantle, we prescribe the same steady, aseismic slip as that in Step 2. On the interface along the top of the subducting oceanic crust and the continental crust and mantle we create two earthquake ruptures. The first rupture applies the coseismic slip from Step 1 at 150 years, while the second rupture prescribes the same steady, aseismic slip as in Step 2.

**pylithapp.problem.formulation.output.domain** Gives the base filename for HDF5 output (for example, `step03.h5`).

We run this example by typing

```{code-block} console
$ pylith step03.cfg
```
The simulation will produce pairs of HDF5/Xdmf files with separate files for each material and fault interface.
{numref}`fig:example:subduction:2d:step03`, which was created using ParaView, displays the magnitude of the displacement field with the deformation exaggerated by a factor of 1000.
Using the animation features within ParaView or Visit you can illustrate how the continental crust near the trench rebounds during the earthquake after subsiding during the interseismic deformation.

:::{figure-md} fig:example:subduction:2d:step03
<img src="figs/step03_soln.*" alt="Solution for Step 3 at 150 years (immediately following the earthquake rupture). The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 3 at 150 years (immediately following the earthquake rupture). The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000.
:::
