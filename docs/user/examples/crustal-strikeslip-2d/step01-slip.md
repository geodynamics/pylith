# Step 1: Static Uniform Coseismic Slip

% Metadata extracted from parameter files.
```{include} step01_slip-synopsis.md
```

## Simulation parameters

This example involves a static simulation that solves for the deformation from prescribed uniform coseismic slip on each of the faults.
{numref}`fig:example:crustal:strikeslip:2d:step01:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step01_slip.cfg`.

:::{figure-md} fig:example:crustal:strikeslip:2d:step01:diagram
<img src="figs/step01-diagram.*" alt="" scale="75%">

Boundary conditions for static coseismic slip.
On the boundary of the domain, we set the displacement component that is perpendicular to the boundary to zero.
We prescribe 4.0 m of right-lateral slip on the main fault, 2.0 m of left-lateral slip on the west branch, and 1.0 m of left-lateral slip on the east branch.
:::

```{code-block} cfg
---
caption: Prescribed slip parameters for Step 1.
---
[pylithapp.problem.interfaces.main_fault.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, -4.0*m, 0.0*m]

[pylithapp.problem.interfaces.west_branch.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, 3.0*m, 0.0*m]

[pylithapp.problem.interfaces.east_branch.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, 1.0*m, 0.0*m]
```

## Running the simulation

```{code-block} console
---
caption: Run Step 1 simulation
---
$ pylith step01_slip.cfg

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
ksp_converged_reason = true
ksp_guess_pod_size = 8
ksp_guess_type = pod
snes_converged_reason = true
snes_monitor = true
ts_error_if_step_fails = true
ts_monitor = true

 >> /src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:239:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const PetscOptions &)
 -- petscoptions(info)
 -- Using user values rather then the following default PETSc options:
ksp_atol = 1.0e-12
ksp_error_if_not_converged = true
ksp_rtol = 1.0e-12
snes_atol = 1.0e-9
snes_error_if_not_converged = true
snes_rtol = 1.0e-12

 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/TimeDependent.py:132:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 4.696586307427e-02
      Linear solve converged due to CONVERGED_ATOL iterations 19
    1 SNES Function norm 7.572811652012e-13
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

At the beginning of the output written to the terminal, we see that PyLith is reading the mesh using the `MeshIOPetsc` reader and that it found the domain to extend from 410000 to 490000 in the x direction and from 3.91e+06 to 3.99e+06 in the y direction.
The scales for nondimensionalization remain the default values for a quasistatic problem.
PyLith detects the presence of a fault based on the Lagrange multiplier for the fault in the solution field and selects appropriate preconditioning options as discussed in {ref}`sec-user-run-pylith-petsc-options`.

At the end of the output written to the terminal, we see that the solver advanced the solution one time step (static simulation).
The linear solve converged after 19 iterations and the norm of the residual met the absolute convergence tolerance (`ksp_atol`) .
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:crustal:strikeslip:2d:step01:solution` we use the `pylith_viz` utility to visualize the y displacement field.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step01_slip-domain.h5 warp_grid --component=y
```

:::{figure-md} fig:example:crustal:strikeslip:2d:step01:solution
<img src="figs/step01-solution.*" alt="Solution for Step 1. The colors indicate the y displacement, and the deformation is exaggerated by a factor of 1000." width="600px"/>

Solution for Step 1.
The colors of the shaded surface indicate the y displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is shown by the gray wireframe.
The contrast in material properties across the faults causes the asymmetry in the y displacement field.
:::

## Step 1 with Cubit Mesh

Using the Cubit mesh rather than the Gmsh mesh involves two changes:

1. Use the `MeshIOCubit` reader instead of the `MeshIOPetsc` reader and change the filename of the mesh file.
2. Set the `label_value` to 1 for boundary conditions and faults.\
   We must override the `label_value` settings in `pylithapp.cfg` that were appropriate for our Gmsh reader but are incorrect for the Cubit reader.

The file `step01_slip_cubit.cfg` provides these changes and updates the names for output.

```{code-block} console
---
caption: Run Step 1 simulation with the Cubit mesh
---
$ pylith step01_slip_cubit.cfg

# The output should look something like the following.
 >> /Users/baagaard/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> /Users/baagaard/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:38:read
 -- meshiocubit(info)
 -- Reading finite-element mesh
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:148:void pylith::meshio::MeshIOCubit::_readVertices(ExodusII &, scalar_array *, int *, int *) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 1610 vertices.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:208:void pylith::meshio::MeshIOCubit::_readCells(ExodusII &, int_array *, int_array *, int *, int *) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 3125 cells in 1 blocks.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:270:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Found 10 node sets.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'boundary_south' with id 10 containing 25 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'boundary_east' with id 11 containing 25 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'boundary_north' with id 12 containing 24 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'boundary_west' with id 13 containing 23 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'fault_main' with id 20 containing 37 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'fault_west' with id 21 containing 13 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'fault_east' with id 22 containing 6 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'fault_main_ends' with id 30 containing 2 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'fault_west_ends' with id 31 containing 2 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:296:void pylith::meshio::MeshIOCubit::_readGroups(ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'fault_east_ends' with id 32 containing 2 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:85:void pylith::meshio::MeshIO::read(pylith::topology::Mesh *, const bool)
 -- meshiocubit(info)
 -- Component 'reader': Domain bounding box:
    (410000, 490000)
    (3.91e+06, 3.99e+06)

-- many lines omitted --

 >> /Users/baagaard/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/TimeDependent.py:132:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 2.540894009864e-02
      Linear solve converged due to CONVERGED_ATOL iterations 45
    1 SNES Function norm 1.131725395538e-12
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The `MeshIOCubit` reader includes diagnostic information in the journal output related to the sizes of the nodesets and material blocks.
