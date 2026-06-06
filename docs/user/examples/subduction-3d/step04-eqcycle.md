# Step 4: Earthquake Cycle with Prescribed Slip

This example combines the interseismic deformation from Step 3 and expands on the earthquake ruptures from Step 2.
We expand the earthquake ruptures to span the entire along strike length of the top of the slab and also consider earthquake rupture on the splay fault.
We use linear Maxwell viscoelastic bulk rheologies in the mantle and deeper part of the slab.
{numref}`fig:example:subduction:3d:step04:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:subduction:3d:step04:diagram
<img src="figs/step04-diagram.*" alt="" width="75%">

Boundary conditions for quasi-static interseismic deformation.
We prescribe aseismic slip (creep) on the bottom of the slab and the deeper portion of the top of the slab; the shallow portion of the top of the slab remains locked.
:::

% Features extracted from simulation parameter files.
```{include} step04_eqcycle-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step04_eqcycle.cfg` and include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for the time step information as well as solution field with displacement and Lagrange multiplier subfields.
* `pylithapp.interfaces` Parameters for the earthquake ruptures and aseismic slip (creep) on the top and bottom of the slab.
* `pylithapp.problem.bc` Parameters for describing the boundary conditions that override the defaults.

We extend the duration of the simulation to 300 years.
We impose two earthquake ruptures on slab interface at t=100 and t=200 years and one earthquake rupture on the splay fault at t=250 years.
Now that we have both earthquake rupture and aseismic creep on the top of the slab, we use `SimpleDB` spatial databases to give depth-dependent, complementary slip.
We impose uniform slip on the splay fault using a `UniformDB`.
As in Step 3 we also adjust the nodesets used for the boundary conditions to remove overlap with the slab to allow the slab to move independently. Note that the aseismic slip uses `KinSrcConstRate`, while the two ruptures use the default `KinSrcStep`.

```{code-block} pyrejournal
---
caption: Run Step 4 simulation using the Gmsh mesh.
---
$ pylith step04_eqcycle.cfg mat_viscoelastic.cfg

# The output should look something like the following.
 >> /software/pylith-opt/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> /software/pylith-opt/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:41:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:74:void pylith::meshio::MeshIO::read(pylith::topology::Mesh*, bool)
 -- meshiopetsc(info)
 -- Component 'meshiopetsc.reader': Domain bounding box:
    (-400000, 400000)
    (-400000, 400000)
    (-400000, 2017.5)

# -- many lines omitted --

 >> /src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:262:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t&, const char*, const pylith::utils::PetscOptions&)
 -- petscoptions(info)
 -- Setting PETSc options:
dm_reorder_section = true
dm_reorder_section_type = cohesive
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_gmres_restart = 100
ksp_guess_pod_size = 8
ksp_guess_type = pod
ksp_rtol = 1.0e-14
mg_fine_ksp_max_it = 5
mg_fine_pc_type = vpbjacobi
pc_type = gamg
snes_atol = 5.0e-7
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-14
ts_error_if_step_fails = true
ts_exact_final_time = matchstep
ts_monitor = true
ts_type = beuler
viewer_hdf5_collective = true

# -- many lines omitted --

30 TS dt 0.1 time 2.9
    0 SNES Function norm 3.515714246561e+03
      Linear solve converged due to CONVERGED_ATOL iterations 22
    1 SNES Function norm 8.157356157583e-07
      Linear solve converged due to CONVERGED_ATOL iterations 1
    2 SNES Function norm 6.740861225135e-07
    Nonlinear solve converged due to CONVERGED_SNORM_RELATIVE iterations 2
31 TS dt 0.1 time 3.
 >> /software/pylith-opt/lib/python3.12/site-packages/pylith/problems/Problem.py:232:finalize
 -- timedependent(info)
 -- Finalizing problem.
 ```

The beginning of the output is near the same as in Steps 2 and 3.
The simulation advances 31 time steps.
As in Step 3 each linear solve requires about 20 iterations to converge.
Some time steps require 2 SNES iterations (likley due to poor mesh quality along the trench) although the linear solve in the second SNES iteration requires no more than 1 iterations in most cases.

```{code-block} pyrejournal
---
caption: Alternatively, run Step 4 simulation using the Cubit mesh.
---
$ pylith step04_eqcycle.cfg step04_eqcycle_cubit.cfg mat_viscoelastic.cfg
```

## Visualizing the results

In {numref}`fig:example:subduction:3d:step04:solution` we use the `pylith_viz` utility to visualize the x displacement field.
You can move the slider or use the `p` and `n` keys to change the increment or decrement time.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step04_eqcycle-domain.h5 warp_grid --component=x --exaggeration=5000
```

:::{figure-md} fig:example:subduction:3d:step04:solution
<img src="figs/step04-solution.*" alt="Solution for Step 4 at t=100 yr. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 5000." width="600px"/>

Solution for Step 4 at t=100 yr.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 5000.
:::
