# Step 2: Earthquake Rupture and Postseismic Relaxation

This example involves a quasi-static simulation for coseismic earthquake rupture and postseismic relaxation.
We use linear Maxwell viscoelastic bulk rheologies in the mantle and deeper part of the slab.
{numref}`fig:example:subduction:3d:step02:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:subduction:3d:step02:diagram
<img src="figs/step02-diagram.*" alt="" width="75%">

Boundary conditions for quasi-static coseismic slip on the subduction interface and postseismic relaxation.
We prescribe uniform oblique slip in the center of the subduction interfaces with roller boundary conditions on the lateral sides and bottom of the domain.
:::

% Features extracted from simulation parameter files.
```{include} step02_coseismic-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step02_coseismic.cfg` and include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for the time step information as well as solution field with displacement and Lagrange multiplier subfields.
* `pylithapp.interfaces` Parameters for the earthquake rupture.

We define the duration of the simulation to be 200 years with an initial time step of 0.1 years.
We enable adaptive time stepping and start with a small time step to resolve the viscoelastic deformation.

For the fault slip, we use the boundary groups for the fault and its buried edges corresponding to the central patch of the top of the slab.
We set the slip initiation time to 10 years with the default step time function for the earthquake rupture.
Alternatively, we could have set the `origin_time` of the earthquake rupture to 10 years and used an `initiation_time` of 0 years.
We impose oblique slip with 1.0 m of right-lateral slip and 4.0 m of reverse slip.

```{code-block} pyrejournal
---
caption: Run Step 2 simulation using the Gmsh mesh.
---
$ pylith step02_coseismic.cfg mat_viscoelastic.cfg

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
ts_adapt_monitor = true
ts_adapt_reject_safety = 0.1
ts_adapt_safety = 0.2
ts_adapt_type = basic
ts_atol = 0.05
ts_error_if_step_fails = true
ts_exact_final_time = matchstep
ts_monitor = true
ts_rtol = 0.05
ts_type = beuler
viewer_hdf5_collective = true

 >> /software/pylith-opt/lib/python3.12/site-packages/pylith/problems/TimeDependent.py:145:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.001 time -0.001
    0 SNES Function norm 2.954073482783e+01
      Linear solve converged due to CONVERGED_ATOL iterations 16
    1 SNES Function norm 1.956521733099e-08
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
      TSAdapt basic beuler 0: step   0 accepted t=-0.001     + 1.000e-03 dt=1.000e-03 

# -- many lines omitted --

9 TS dt 0.496292 time 1.50371
    0 SNES Function norm 2.765617769230e+00
      Linear solve converged due to CONVERGED_ATOL iterations 12
    1 SNES Function norm 1.288354163470e-07
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
      TSAdapt basic beuler 0: step   9 accepted t=1.50371    + 4.963e-01 dt=1.373e+00  wlte=0.00523  wltea=   -1 wlter=   -1
9 TS dt 0. time 1.50371
    0 SNES Function norm 2.765639325307e+00
      Linear solve converged due to CONVERGED_ATOL iterations 14
    1 SNES Function norm 9.313028842369e-09
    Nonline496294ar solve converged due to CONVERGED_FNORM_ABS iterations 1
      TSAdapt basic beuler 0: step   9 accepted t=1.50371    + 4.963e-01 dt=1.373e+00  wlte=0.00523  wltea=   -1 wlter=   -1
10 TS dt 1.37262 time 2.
 >> /software/pylith-opt/lib/python3.12/site-packages/pylith/problems/Problem.py:232:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

At the beginning of the output written to the terminal, we see that PyLith is reading the mesh using the `MeshIOPetsc` reader.
We also see the PETSc solver options, which show use of the variable point-block Jacobi and GAMG (algebriac multigrid) preconditioner settings.

At the end of the output written to the terminal, we see that the solver advanced the solution 10 time steps.
The linear solve converged after 14 iterations and the norm of the residual met the absolute convergence tolerance (`ksp_atol`) .
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).

```{code-block} pyrejournal
---
caption: Alternatively, run Step 2 simulation using the Cubit mesh.
---
$ pylith step02_coseismic.cfg step02_coseismic_cubit.cfg mat_viscoelastic.cfg
```

## Visualizing the results

In {numref}`fig:example:subduction:3d:step02:solution` we use the `pylith_viz` utility to visualize the x displacement field.
You can move the slider or use the `p` and `n` keys to change the increment or decrement time.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step02_coseismic-domain.h5 warp_grid --component=x --exaggeration=10000
```

:::{figure-md} fig:example:subduction:3d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2 at t=200 yr. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 10,000." width="600px"/>

Solution for Step 2 at t=200 yr.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 10,000.
:::
