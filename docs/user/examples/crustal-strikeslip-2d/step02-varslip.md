# Step 2: Static Spatially Variable Coseismic Slip

% Metadata extracted from parameter files.
```{include} step02_varslip-synopsis.md
```

## Simulation parameters

This example involves a static simulation that solves for the deformation from prescribed spatially variable coseismic slip on each of the faults.
{numref}`fig:example:crustal:strikeslip:2d:step02:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step02_varslip.cfg`.

:::{figure-md} fig:example:crustal:strikeslip:2d:step02:diagram
<img src="figs/step02-diagram.*" alt="" scale="75%">

Boundary conditions for static coseismic slip.
On the boundary of the domain, we set the displacement component that is perpendicular to the boundary to zero.
We prescribe piecewise linear variations in slip along strike on each of the faults.
:::

We prescribe piecewise linear variations in slip along strike on each of the faults using `SimpleDB` spatial databases.
The spatial databases contain coordinates of points in geographic coordinates.
PyLith will transform each the coordinates of point on the fault to the geographic coordinates and interpolate the values in the spatial database to the point on the fault.

```{code-block} cfg
---
caption: Prescribed slip parameters for Step 2.
---
[pylithapp.problem.interfaces.main_fault.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Slip parameters for fault 'main'
db_auxiliary_field.iohandler.filename = fault_main_slip.spatialdb
db_auxiliary_field.query_type = linear

[pylithapp.problem.interfaces.west_branch.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Slip parameters for fault 'west'
db_auxiliary_field.iohandler.filename = fault_west_slip.spatialdb
db_auxiliary_field.query_type = linear

[pylithapp.problem.interfaces.east_branch.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Slip parameters for fault 'east'
db_auxiliary_field.iohandler.filename = fault_east_slip.spatialdb
db_auxiliary_field.query_type = linear
```

## Running the simulation

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02_varslip.cfg

# The output should look something like the following.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:38:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:85:void pylith::meshio::MeshIO::read(pylith::topology::Mesh *, const bool)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (410000, 490000)
    (3.91e+06, 3.99e+06)
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'fault_main'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'fault_west'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'fault_east'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:116:preinitialize
 -- timedependent(info)
 -- Performing minimal initialization before verifying configuration.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Solution.py:39:preinitialize
 -- solution(info)
 -- Performing minimal initialization of solution.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'fault_main'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'fault_west'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/faults/FaultCohesiveKin.py:87:preinitialize
 -- faultcohesivekin(info)
 -- Pre-initializing fault 'fault_east'.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:174:verifyConfiguration
 -- timedependent(info)
 -- Verifying compatibility of problem configuration.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:219:_printInfo
 -- timedependent(info)
 -- Scales for nondimensionalization:
    Length scale: 1000*m
    Time scale: 3.15576e+09*s
    Pressure scale: 3e+10*m**-1*kg*s**-2
    Density scale: 2.98765e+23*m**-3*kg
    Temperature scale: 1*K
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:185:initialize
 -- timedependent(info)
 -- Initializing timedependent problem with quasistatic formulation.
 >> /src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:239:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const PetscOptions &)
 -- petscoptions(info)
 -- Setting PETSc options:
dm_reorder_section = true
dm_reorder_section_type = cohesive
ksp_atol = 1.0e-7
ksp_converged_reason = true
ksp_error_if_not_converged = true
ksp_guess_pod_size = 8
ksp_guess_type = pod
ksp_rtol = 1.0e-12
mg_fine_ksp_max_it = 5
mg_fine_pc_type = vpbjacobi
mg_levels_pc_type = pbjacobi
pc_gamg_coarse_eq_limit = 200
pc_type = gamg
snes_atol = 4.0e-7
snes_converged_reason = true
snes_error_if_not_converged = true
snes_monitor = true
snes_rtol = 1.0e-12
ts_error_if_step_fails = true
ts_exact_final_time = matchstep
ts_monitor = true
ts_type = beuler
viewer_hdf5_collective = true

 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/TimeDependent.py:132:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.001 time 0.
    0 SNES Function norm 1.139546193525e-01
      Linear solve converged due to CONVERGED_ATOL iterations 14
    1 SNES Function norm 1.249276796368e-09
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.001 time 0.001
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

At the end of the output written to the terminal, we see that the solver advanced the solution one time step (static simulation).
The linear solve converged after 14 iterations and the norm of the residual met the absolute convergence tolerance (`ksp_atol`) .
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).

## Visualizing the results

In {numref}`fig:example:crustal:strikeslip:2d:step02:solution` we use the `pylith_viz` utility to visualize the y displacement field.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step02_varslip-domain.h5 warp_grid --component=y
```

:::{figure-md} fig:example:crustal:strikeslip:2d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2. The colors indicate the y displacement, and the deformation is exaggerated by a factor of 1000." width="600px"/>

Solution for Step 2.
The colors of the shaded surface indicate the y displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is shown by the gray wireframe.
The contrast in material properties across the faults causes the asymmetry in the y displacement field.
:::
