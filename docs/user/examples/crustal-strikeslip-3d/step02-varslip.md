# Step 2: Static Spatially Variable Coseismic Slip

% Metadata extracted from parameter files.
```{include} step02_varslip-synopsis.md
```

## Simulation parameters

This example involves a static simulation that solves for the deformation from prescribed spatially varying coseismic slip on each of the faults.
{numref}`fig:example:crustal:strikeslip:3d:step02:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step02_varslip.cfg`.

:::{figure-md} fig:example:crustal:strikeslip:3d:step02:diagram
<img src="figs/step02-diagram.*" alt="" scale="50%">

Boundary conditions for static coseismic slip.
On the boundary of the domain, we set the displacement component that is perpendicular to the boundary to zero.
We prescribe slip that varies linearly over each of the fault surfaces.
:::

We prescribe piecewise linear variations in slip along strike on each of the faults using `SimpleDB` spatial databases
The spatial databases contain coordinates of points in geographic coordinates.
PyLith will transform each the coordinates of point on the fault to the geographic coordinates and interpolate the values in the spatial database to the point on the fault.

:::{tip}
In most cases you will not want the interpolation to be done in geographic coordinates, because distances in each of the coordinate directions are different.
This is especially true when horizontal coordinates are given in longitude and latitude and vertical distances are given in meters.
Transforming the points to a geographic projection is generally preferred.
:::

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

```{code-block} pyrejournal
---
caption: Run Step 2 simulation
---
$ pylith step02_varslip.cfg

 >> software/pylith-opt/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:76:main
 -- info (application-flow)
 -- Running on 1 process(es).
 >> src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:251:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t&, const char*, const pylith::utils::PetscOptions&)
 -- info (application-flow)
 -- Setting PETSc options:
    dm_reorder_section = true
    dm_reorder_section_type = cohesive
    ksp_atol = 1.0e-7
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

 >> src/cig/pylith/libsrc/pylith/meshio/MeshIOPetsc.cc:204:virtual void pylith::meshio::MeshIOPetsc::_read()
 -- info (application-flow)
 -- Component 'meshiopetsc.reader': Reading finite-element mesh from 'mesh_tet.msh'.
 >> src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:76:void pylith::meshio::MeshIO::read(pylith::topology::Mesh*, bool)
 -- info (application-flow)
 -- Component 'meshiopetsc.reader': Domain bounding box:
    (413700, 493700)
    (3.917e+06, 3.977e+06)
    (-40000, 0)
 >> src/cig/pylith/libsrc/pylith/initializers/MeshDistributor.cc:57:virtual pylith::topology::Mesh* pylith::initializers::MeshDistributor::run(pylith::topology::Mesh*, const pylith::problems::Problem&)
 -- info (application-flow)
 -- Distributing mesh.
 >> src/cig/pylith/libsrc/pylith/initializers/MeshInsertInterfaces.cc:50:virtual pylith::topology::Mesh* pylith::initializers::MeshInsertInterfaces::run(pylith::topology::Mesh*, const pylith::problems::Problem&)
 -- info (application-flow)
 -- Inserting cohesive cells.
 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:316:virtual void pylith::problems::TimeDependent::verifyConfiguration() const
 -- info (application-flow)
 -- Component 'timedependent.problem': Verifying problem configuration.
 >> software/pylith-opt/lib/python3.12/site-packages/pylith/problems/Problem.py:238:_printInfo
 -- info (application-flow)
 -- Scales for nondimensionalization:
    Length scale: 2500*m
    Displacement scale: 1*m
    Time scale: 3.15576e+09*s
    Rigidity scale: 1e+10*m**-1*kg*s**-2
    Temperature scale: 1*K
 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:342:virtual void pylith::problems::TimeDependent::initialize()
 -- info (application-flow)
 -- Component 'timedependent.problem': Initializing problem.
 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:473:void pylith::problems::TimeDependent::solve()
 -- info (application-flow)
 -- Component 'timedependent.problem': Solving equations.
0 TS dt 0.001 time 0.
    0 SNES Function norm 1.091558353128e+01
      Linear solve converged due to CONVERGED_ATOL iterations 18
    1 SNES Function norm 2.786048497788e-08
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.001 time 0.001
 >> software/pylith-opt/lib/python3.12/site-packages/pylith/problems/Problem.py:222:finalize
 -- info (application-flow)
 -- Finalizing problem.
```

At the beginning of the output written to the terminal, we see that PyLith is reading the mesh using the `MeshIOPetsc` reader and that it found the domain to extend from 413700 to 493700 in the x direction, from 3.92e+06 to 3.98e+06 in the y direction, and from -40000 to 0 in the z direction.
The scales for nondimensionalization remain the default values for a quasistatic problem.
PyLith detects the presence of a fault based on the Lagrange multiplier for the fault in the solution field and selects appropriate preconditioning options as discussed in {ref}`sec-user-run-pylith-petsc-options`.

At the end of the output written to the terminal, we see that the solver advanced the solution one time step (static simulation).
The linear solve converged after 18 iterations and the norm of the residual met the absolute convergence tolerance (`ksp_atol`) .
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).

## Visualizing the results

In {numref}`fig:example:crustal:strikeslip:3d:step02:solution` we use the `pylith_viz` utility to visualize the y displacement field.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step02_varslip-domain.h5 warp_grid --component=y
```

:::{figure-md} fig:example:crustal:strikeslip:3d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2. The colors indicate the y displacement, and the deformation is exaggerated by a factor of 1000." width="600px"/>

Solution for Step 2.
The colors of the shaded surface indicate the y displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is shown by the gray wireframe.
The contrast in material properties across the faults causes the asymmetry in the y displacement field.
:::
