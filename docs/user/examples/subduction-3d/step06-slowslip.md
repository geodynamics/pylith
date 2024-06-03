# Step 6: Prescribed Slow Slip Events

This example simulates several slow slip events on the slab interface.
We prescribe a time history of slip, varying the slip amplitude as a function of time.
We also output the displacement field at fake CGNSS stations using `OutputSolnPoints`.
For this simulation, we use linearly elastic materials.
{numref}`fig:example:subduction:3d:step06:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:subduction:3d:step06:diagram
<img src="figs/step06-diagram.*" alt="" width="75%">

Boundary conditions for prescribed slow slip.
We prescribe time-varying slip on a patch within the subduction interface with roller boundary conditions on the lateral sides and bottom of the domain.
:::

% Features extracted from simulation parameter files.
```{include} step06_slowslip-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step06_slowslip.cfg` and include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for the time step information as well as solution field with displacement and Lagrange multiplier subfields.
* `pylithapp.interfaces` Parameters for the time-varying slip on the top of the slab.

For slow slip we use the `KinSrcTimeHistory` kinematic source to prescribe time-varying slip. We also use `OutputSolnPoints` to simulate output at CGNSS stations.

Prior to running the simulation, we use the `generate_slowslip.py` Python script, which in turn reads the parameters in `generate_slowslip.cfg`. This will generate a spatial database and a time database needed for the simulation:
```{code-block} console
---
caption: Generate database files needed for Step 6
---
# Generate fault_slabtop_slowslip.spatialdb and fault_slabtop_slowslip.timedb.
$ ./generate_slowslip.py
```
Once the database files have been generated we can run the simulation.

```{code-block} console
---
caption: Run Step 6 simulation
---
$ pylith step06_slowslip.cfg mat_elastic.cfg

# The output should look something like the following.
 >> software/virtualenv/python312/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> software/virtualenv/python312/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:38:read
 -- meshiocubit(info)
 -- Reading finite-element mesh
 >> software/cig/pylith3/source/pylith-fork/libsrc/pylith/meshio/MeshIOCubit.cc:148:void pylith::meshio::MeshIOCubit::_readVertices(ExodusII &, scalar_array *, int *, int *) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 24738 vertices.
 >> software/cig/pylith3/source/pylith-fork/libsrc/pylith/meshio/MeshIOCubit.cc:208:void pylith::meshio::MeshIOCubit::_readCells(ExodusII &, int_array *, int_array *, int *, int *) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 133827 cells in 4 blocks.

# -- many lines omitted --

 >> software/cig/pylith3/source/pylith-fork/libsrc/pylith/utils/PetscOptions.cc:239:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const PetscOptions &)
 -- petscoptions(info)
 -- Setting PETSc options:
dm_reorder_section = true
dm_reorder_section_type = cohesive
ksp_atol = 1.0e-12
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_guess_pod_size = 8
ksp_guess_type = pod
ksp_rtol = 1.0e-12
mg_fine_pc_type = vpbjacobi
pc_type = gamg
snes_atol = 1.0e-9
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_monitor = true
ts_type = beuler

# -- many lines omitted --

15 TS dt 5.4757e-05 time 0.000766598
    0 SNES Function norm 2.850925995479e-01
      Linear solve converged due to CONVERGED_ATOL iterations 4
    1 SNES Function norm 2.450595224581e-10
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
16 TS dt 5.4757e-05 time 0.000821355
    0 SNES Function norm 2.450592912419e-10
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 0
17 TS dt 5.4757e-05 time 0.000876112
 >> software/virtualenv/python312/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.

```

The beginning of the output is nearly the same as in previous examples. The simulation advances 17 time steps; however, no slip occurs after step 15.

## Visualizing the results

In {numref}`fig:example:subduction:3d:step06:solution` we use the `pylith_viz` utility to visualize the x displacement field.
You can move the slider or use the `p` and `n` keys to change the increment or decrement time.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step06_slowslip-domain.h5 warp_grid --component=x --exaggeration=10000
```

:::{figure-md} fig:example:subduction:3d:step06:solution
<img src="figs/step06-solution.*" alt="Solution for Step 6 at t=0.09 yr. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 10,000." width="600px"/>

Solution for Step 6 at t=0.09 yr.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 10,000.
:::
