# Step 8: Stress Field Due to Gravitational Body Forces

This example demonstrates the use of gravitational body forces as well as the use of initial stresses to balance the body forces.
This involves enabling gravity within our domain with Dirichlet roller boundary conditions on the lateral and bottom boundaries; we do not include faults in this example.
We also demonstrate what happens when the initial stresses are not in balance with the gravitational stresses, and show how viscoelastic problems with gravitational stresses will in general not reach a steady-state solution.
The example is divided into three sub-problems:

**Step 8a**  Gravitational body forces with 3-D density variations in elastic materials and initial stresses for a uniform density.

**Step 8b**  Gravitational body forces with 3-D density variations in elastic materials and initial stresses from Step 8a (initial stresses satisfy equilibrium, so there is almost no deformation).

**Step 8c**  Gravitational body forces with 3-D density variations in elastic and viscoelastic materials and initial stresses from Step 8a plus finite strain formulation (does not reach a steady-state solution).

## Step 08a

For Step 8a we apply gravitational stresses and attempt to balance these with analytically computed stresses consistent with the density of the mantle.
Since the actual density is not uniform, the stresses are out of balance and we end up with some deformation.
In `step08a.cfg` we turn on gravity and set the total time to zero (there is no time dependence in this model).

```{code-block} cfg
---
caption: Excerpt from `step08a.cfg`
---
[pylithapp.problem]
# Set gravity field (default is None).
gravity_field = spatialdata.spatialdb.GravityField

[pylithapp.problem.formulation.time_step]
# Define the total time for the simulation.
total_time = 0.0*year
```

Our initial stress field corresponds to {math}`\sigma_{xx} = \sigma_{yy} = \sigma_{zz} = \rho_{mantle}gz` for all four materials, where {math}`mantle` is the density of the mantle, {math}`g` is the acceleration due to gravity, and {math}`z` is elevation.
With only two control points necessary to describe this linear variation, we use the same *SimpleDB* spatial database for all four materials.

```{code-block} cfg
---
caption: Excerpt from `step08a.cfg`
---
# We specify initial stresses for each material via a SimpleDB and linear interpolation.
[pylithapp.problem.materials.slab]
db_initial_stress = spatialdata.spatialdb.SimpleDB
db_initial_stress.label = Initial stress in the slab
db_initial_stress.iohandler.filename = spatialdb/mat_initial_stress_grav.spatialdb
db_initial_stress.query_type = linear

[pylithapp.problem.materials.wedge]
db_initial_stress = spatialdata.spatialdb.SimpleDB
db_initial_stress.label = Initial stress in the wedge
db_initial_stress.iohandler.filename = spatialdb/mat_initial_stress_grav.spatialdb
db_initial_stress.query_type = linear

[pylithapp.problem.materials.mantle]
db_initial_stress = spatialdata.spatialdb.SimpleDB
db_initial_stress.label = Initial stress in the mantle
db_initial_stress.iohandler.filename = spatialdb/mat_initial_stress_grav.spatialdb
db_initial_stress.query_type = linear

[pylithapp.problem.materials.crust]
db_initial_stress = spatialdata.spatialdb.SimpleDB
db_initial_stress.label = Initial stress in the crust
db_initial_stress.iohandler.filename = spatialdb/mat_initial_stress_grav.spatialdb
db_initial_stress.query_type = linear
```

```{code-block} console
---
caption: Run Step 8a simulation
---
$ pylith step08a.cfg mat_elastic.cfg solver_algebraicmultigrid.cfg
```

The simulation will generate ten pairs of HDF5/Xdmf files beginning with `step08a`:

**step08a-domain.h5[.xmf]**  Time series of the solution field over the domain.

**step08a-groundsurf.h5[.xmf]**  Time series of the solution field over the ground surface.

**step08a-slab_info.h5[.xmf]**  Properties for the slab material.

**step08a-slab.h5[.xmf]**  Time series of the state variables (stress and strain) for the slab material.

**step08a-wedge_info.h5[.xmf]**  Properties for the wedge material.

**step08a-wedge.h5[.xmf]**  Time series of the state variables (stress and strain) for the wedge material.

**step08a-crust_info.h5[.xmf]**  Properties for the crust material.

**step08a-crust.h5[.xmf]**  Time series of the tate variables (stress and strain) for the crust material.

**step08a-mantle_info.h5[.xmf]**  Properties for the mantle material.

**step08a-mantle.h5[.xmf]**  Time series of the state variables (stress and strain) for the mantle material.

When the problem has run, we see deformation that is consistent with the mismatched densities.
The slab subsides while the crust undergoes uplift due to the differences in density relative to the mantle.
{numref}`fig:example:subduction:3d:step08a` shows the deformed mesh visualized with the `plot_dispwarp.py` ParaView Python script.

:::{figure-md} fig:example:subduction:3d:step08a
<img src="figs/subduction3d_step08a_soln.*" alt="Solution for Step 8a. The deformation has been exaggerated by a factor of 500 and the colors highlight the vertical displacement component. The crustal material in the east is less dense than the assumed mantle material for initial stresses, while the slab material in the west is more dense. The result is uplift in the east and subsidence in the west." width="100%"/>

Solution for Step 8a. The deformation has been exaggerated by a factor of 500 and the colors highlight the vertical displacement component. The crustal material in the east is less dense than the assumed mantle material for initial stresses, while the slab material in the west is more dense. The result is uplift in the east and subsidence in the west.
:::

## Step 8b

Step 8b is similar to Step 8a, but we use the stresses output from Step 8a as the initial stress rather than analytically computing initial stresses.
Because the initial stresses are consistent with the variations in density across the materials, the initial stresses will satisfy equilibrium and there will be essentially no deformation; the initial stresses do not perfectly balance because in Step 8a we average the values over the quadrature points for the output.
We use the Python script `generate_initial_stress.py`, located in the `spatialdb` directory, to postprocess the output from Step 8a and generate the initial stress spatial database.
Note that this script uses the Python interface to the spatialdata package to write the spatial database; this is much easier than writing a script to format the data to conform to the format of the spatial database.
The spatial database will contain the stresses at each cell of our unstructured mesh, so the points are not on a logical grid, and we must use a *SimpleDB*.

```{code-block} console
---
caption: Generate the initial stresses for Step 8b
---
# From the examples/3d/subduction directory, change to the spatialdb subdirectory.
$ cd spatialdb
$ ./generate_initial_stress.py
```

This will create spatial databases containing initial stresses for each of the four materials.

In the `step08b.cfg` file we specify the *SimpleDB* spatial database for each material (they are now material specific).
With points at each cell centroid, we use nearest interpolation (default) rather than linear interpolation; this is a small approximation but it is much faster than using linear interpolation in this unstructured set of points.

```{code-block} cfg
---
caption: Excerpt from `step08b.cfg`
---
[pylithapp.problem.materials.slab]
db_initial_stress = spatialdata.spatialdb.SimpleDB
db_initial_stress.label = Initial stress in the slab
db_initial_stress.iohandler.filename = spatialdb/mat_initial_stress_grav-slab.spatialdb

[pylithapp.problem.materials.wedge]
db_initial_stress = spatialdata.spatialdb.SimpleDB
db_initial_stress.label = Initial stress in the wedge
db_initial_stress.iohandler.filename = spatialdb/mat_initial_stress_grav-wedge.spatialdb

[pylithapp.problem.materials.mantle]
db_initial_stress = spatialdata.spatialdb.SimpleDB
db_initial_stress.label = Initial stress in the mantle
db_initial_stress.iohandler.filename = spatialdb/mat_initial_stress_grav-mantle.spatialdb

[pylithapp.problem.materials.crust]
db_initial_stress = spatialdata.spatialdb.SimpleDB
db_initial_stress.label = Initial stress in the crust
db_initial_stress.iohandler.filename = spatialdb/mat_initial_stress_grav-crust.spatialdb
```

```{code-block} console
---
caption: Run Step 8b simulation
---
$ pylith step08b.cfg mat_elastic.cfg solver_algebraicmultigrid.cfg
```
This simulation will produce files in the `output` directory analogous to Step 8a.

When we compare the resulting elastic displacements with those of Step 8a, we find that there is essentially no displacement,
as seen in {numref}`fig:example:subduction:3d:step08b`.

:::{figure-md} fig:example:subduction:3d:step08b
<img src="figs/subduction3d_step08b_soln.*" alt="Solution for Step 8b. In this case the initial stresses satisfy the governing equation, so there is no deformation." width="100%"/>

Solution for Step 8b. In this case the initial stresses satisfy the governing equation, so there is no deformation.
:::

## Step 8c

In this example we use linear Maxwell viscoelastic models in place of the elastic models for the slab and mantle.
We also use the small strain formulation (*ImplicitLgDeform*) so that the deformed configuration is taken into account; Steps 8a and 8b use the default *Implicit* infinitesimal strain formulation.
The small strain formulation should generally be used for viscoelastic problems with gravity where you need accurate estimates of the vertical deformation.

:::{warning}
The shear stress variations in the initial stresses will cause the viscoelastic materials to drive viscous flow, resulting in time-dependent deformation.
As long as the elastic materials impose deviatoric stresses in the viscoelastic materials through continuity of strain, the viscoelastic materials will continue to flow.
**As a result, in this case and many other simulations with viscoelastic materials and gravitational body forces, it is difficult to find a steady state solution.**
:::

The only difference between the parameters in `step08b.cfg` and `step08c.cfg` is in the formulation setting and the simulation time:

```{code-block} cfg
---
caption: Excerpt from `step08c.cfg`
---
[pylithapp.timedependent]
# Turn on the small strain formulation, which automatically runs the
# simulation as a nonlinear problem.
formulation = pylith.problems.ImplicitLgDeform

# Set gravity field (default is None).
gravity_field = spatialdata.spatialdb.GravityField

[pylithapp.problem.formulation.time_step]
# Define the total time for the simulation and the time step size.
total_time = 100.0*year
dt = 10.0*year
```

We use the material settings in `mat_viscoelastic.cfg`.

```{code-block} console
---
caption: Run Step 8c simulation
---
pylith step08c.cfg mat_viscoelastic.cfg
solver_algebraicmultigrid.cfg
```

This simulation will produce files in the `output` directory analogous to Steps 8a and 8b.

The resulting deformation is shown in {numref}`fig:example:subduction:3d:step08c`.
As a result of viscous flow, the vertical deformation is even larger than that for Step 8a.
If we were to run the simulation for a longer time period, the amount of vertical deformation would continue to increase.

:::{figure-md} fig:example:subduction:3d:step08c
<img src="figs/subduction3d_step08c_soln.*" alt="Image generated by running the `plot_dispwarp.py` script for sub-problem step08c. Although the stresses balance in the elastic solution, viscous flow in subsequent time steps results in large vertical deformation." width="100%"/>

Image generated by running the `plot_dispwarp.py` script for sub-problem step08c. Although the stresses balance in the elastic solution, viscous flow in subsequent time steps results in large vertical deformation.
:::

## Exercises

* What happens in sub-problem step08a if we use a different reference density to compute our initial stresses?
* For sub-problem step08b, what happens if, for one of the materials you use the initial stresses from sub-problem step08a?
* For sub-problem step08c, what happens if you:
  * Run the simulation for a longer period of time?
  * Change the viscoelastic properties? For example, reduce the viscosity, make all materials viscoelastic, switch to a power-law rheology, etc.
* Is it possible to find a better initial stress state for sub-problem step08c?
  * What if the initial stresses were computed with nearly incompressible materials, and all materials in the model are viscoelastic?
