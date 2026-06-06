# Step 3: Interseismic Deformation

This example involves a quasi-static simulation for interseismic deformation.
We prescribe aseismic slip (creep) on the bottom of the slab and the deeper portion of the top of the slab; the shallow portion of the top of the slab remains locked.
We use linear Maxwell viscoelastic bulk rheologies in the mantle and deeper part of the slab.
{numref}`fig:example:subduction:3d:step03:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:subduction:3d:step03:diagram
<img src="figs/step03-diagram.*" alt="" width="75%">

Boundary conditions for quasi-static interseismic deformation.
We prescribe aseismic slip (creep) on the bottom of the slab and the deeper portion of the top of the slab; the shallow portion of the top of the slab remains locked.
:::

% Features extracted from simulation parameter files.
```{include} step03_interseismic-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step03_interseismic.cfg` and include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for the time step information as well as solution field with displacement and Lagrange multiplier subfields.
* `pylithapp.interfaces` Parameters for the aseismic slip (creep) on the top and bottom of the slab.
* `pylithapp.problem.bc` Parameters for describing the boundary conditions that override the defaults.

For aseismic slip we use the `KinSrcConstRate` kinematic source to prescribe a constant slip rate.
We also adjust the groups used for the boundary conditions to remove overlap with the slab to allow the slab to move independently.
We adjust the displacement scale to match the slow deformation and tighten the KSP absolute tolerance to facilitate convergence in a single SNES iteration.

```{code-block} pyrejournal
---
caption: Run Step 3 simulation using the Gmsh mesh.
---
$ pylith step03_interseismic.cfg mat_viscoelastic.cfg

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

20 TS dt 0.1 time 1.9
    0 SNES Function norm 3.124244823506e+03
      Linear solve converged due to CONVERGED_ATOL iterations 22
    1 SNES Function norm 1.157676866316e-07
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
21 TS dt 0.1 time 2.
 >> /software/pylith-opt/lib/python3.12/site-packages/pylith/problems/Problem.py:232:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The beginning of the output is nearly the same as in Step 2.
The simulation advances 21 time steps.
The linear solve converged after 22 iterations and the norm of the residual met the absolute convergence tolerance (`ksp_atol`) .
In this simulation the fault interfaces on the top and bottom of the slab occupy a significant fraction of the domain.
As a result, the linear solver requires more iterations to converge compared to the limited fault interface in Step 2.


```{code-block} pyrejournal
---
caption: Alternatively, run Step 3 simulation using the Cubit mesh.
---
$ pylith step03_interseismic.cfg step03_interseismic_cubit.cfg mat_viscoelastic.cfg
```

## Visualizing the results

In {numref}`fig:example:subduction:3d:step03:solution` we use the `pylith_viz` utility to visualize the x displacement field.
You can move the slider or use the `p` and `n` keys to change the increment or decrement time.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step03_interseismic-domain.h5 warp_grid --component=x --exaggeration=5000
```

:::{figure-md} fig:example:subduction:3d:step03:solution
<img src="figs/step03-solution.*" alt="Solution for Step 3 at t=140 yr. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 5000." width="600px"/>

Solution for Step 3 at t=140 yr.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 5000.
:::
