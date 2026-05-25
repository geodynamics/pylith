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
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:79:main
 -- info (application-flow)
 -- Running on 1 process(es).
 >> src/cig/pylith/libsrc/pylith/utils/PetscOptions.cc:251:static void pylith::utils::_PetscOptions::write(pythia::journal::info_t &, const char *, const PetscOptions &)
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

 >> src/cig/pylith/libsrc/pylith/meshio/MeshIOPetsc.cc:205:virtual void pylith::meshio::MeshIOPetsc::_read()
 -- info (application-flow)
 -- Component 'meshiopetsc.reader': Reading finite-element mesh from 'mesh_tri.msh'.
 >> src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:76:void pylith::meshio::MeshIO::read(pylith::topology::Mesh *, const bool)
 -- info (application-flow)
 -- Component 'meshiopetsc.reader': Domain bounding box:
    (410000, 490000)
    (3.91e+06, 3.99e+06)
 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:316:virtual void pylith::problems::TimeDependent::verifyConfiguration() const
 -- info (application-flow)
 -- Component 'timedependent.problem': Verifying problem configuration.
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:238:_printInfo
 -- info (application-flow)
 -- Scales for nondimensionalization:
    Length scale: 1000*m
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
    0 SNES Function norm 2.419639130581e+01
      Linear solve converged due to CONVERGED_ATOL iterations 13
    1 SNES Function norm 6.875321480887e-08
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.001 time 0.001
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:222:finalize
 -- info (application-flow)
 -- Finalizing problem.
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/utils/PetscManager.py:60:finalize
 -- info (application-flow)
 -- Finalizing PETSc.
```

At the beginning of the output written to the terminal, we see that PyLith is reading the mesh using the `MeshIOPetsc` reader and that it found the domain to extend from 410000 to 490000 in the x direction and from 3.91e+06 to 3.99e+06 in the y direction.
PyLith detects the presence of a fault based on the Lagrange multiplier for the fault in the solution field and selects appropriate preconditioning options as discussed in {ref}`sec-user-run-pylith-petsc-options`.
The scales for nondimensionalization remain the default values for a quasistatic problem.

At the end of the output written to the terminal, we see that the solver advanced the solution one time step (static simulation).
The linear solve converged after 13 iterations and the norm of the residual met the absolute convergence tolerance (`ksp_atol`) .
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).

### Earthquake rupture parameters

We use the `pylith_eqinfo` utility to compute rupture information, such as earthquake magnitude, seismic moment, seismic potency, and average slip.
For 2D simulations, the average slip provides the most useful information.
[`pylith_eqinfo`](../../run-pylith/utilities.md) is a Pyre application and you specify parameters using `cfg` files and the command line.
It writes results to a Python script for use in post-processing of simulation results.

The file `eqinfoapp.cfg` holds the parameters for `pylith_eqinfo` in this example.
By default, `pylith_eqinfo` extracts information for the final time step in the output file.
In order to compute the seismic moment and moment magnitude, we need the shear modulus, which in this case is uniform over the domain.
Consquently, we specify the density and shear wave speed using a `UniformDB` in `eqinfoapp.cfg`.

```{code-block} console
---
caption: Run `pylith_eqinfo` for the Step 1 simulation.
---
$ pylith_eqinfo
```

```{code-block} python
---
caption: Contents of `output/step01_slip-eqinfo.py` generated by `pylith_eqinfo`. The `all` object includes earthquake rupture information combined from the three the individual faults. The average slip matches the uniform slip prescribed on each fault.
---
class RuptureStats(object):
    pass
all = RuptureStats()
all.timestamp = [  3.155760e+07]
all.ruparea = [  5.387992e+04]
all.potency = [  1.733080e+05]
all.moment = [  3.899430e+15]
all.avgslip = [  3.216560e+00]
all.mommag = [  4.360667e+00]
main_fault = RuptureStats()
main_fault.timestamp = [  3.155760e+07]
main_fault.ruparea = [  3.577375e+04]
main_fault.potency = [  1.430950e+05]
main_fault.moment = [  3.219637e+15]
main_fault.avgslip = [  4.000000e+00]
main_fault.mommag = [  4.305205e+00]
east_branch = RuptureStats()
east_branch.timestamp = [  3.155760e+07]
east_branch.ruparea = [  5.999357e+03]
east_branch.potency = [  5.999357e+03]
east_branch.moment = [  1.349855e+14]
east_branch.avgslip = [  1.000000e+00]
east_branch.mommag = [  3.386858e+00]
west_branch = RuptureStats()
west_branch.timestamp = [  3.155760e+07]
west_branch.ruparea = [  1.210682e+04]
west_branch.potency = [  2.421364e+04]
west_branch.moment = [  5.448068e+14]
west_branch.avgslip = [  2.000000e+00]
west_branch.mommag = [  3.790828e+00]
```

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
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:79:main
 -- info (application-flow)
 -- Running on 1 process(es).
 
 # -- many lines omitted --

 >> src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:170:virtual void pylith::meshio::MeshIOCubit::_read()
 -- info (application-flow)
 -- Component 'meshiocubit.reader': Reading finite-element mesh from 'mesh_tri.exo'.
 >> src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:207:virtual void pylith::meshio::MeshIOCubit::_read()
 -- info (application-flow)
 -- Component 'meshiocubit.reader': Read 3125 cells and 1610 vertices.

# -- many lines omitted --

>> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:473:void pylith::problems::TimeDependent::solve()
 -- info (application-flow)
 -- Component 'timedependent.problem': Solving equations.
0 TS dt 0.001 time 0.
    0 SNES Function norm 2.428257473425e+01
      Linear solve converged due to CONVERGED_ATOL iterations 14
    1 SNES Function norm 1.922622064928e-08
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.001 time 0.001
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:222:finalize
 -- info (application-flow)
 -- Finalizing problem.
```
