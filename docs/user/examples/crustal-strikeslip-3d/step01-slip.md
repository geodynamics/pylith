# Step 1: Static Coseismic Slip

% Metadata extracted from parameter files.
```{include} step01_slip-synopsis.md
```

## Simulation parameters

This example involves a static simulation that solves for the deformation from prescribed uniform coseismic slip on each of the faults.
{numref}`fig:example:crustal:strikeslip:3d:step01:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step01_slip.cfg`.

:::{figure-md} fig:example:crustal:strikeslip:3d:step01:diagram
<img src="figs/step01-diagram.*" alt="" scale="50%">

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
    (413700, 493700)
    (3.917e+06, 3.977e+06)
    (-40000, 0)

-- many lines omitted --

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
    0 SNES Function norm 1.964459680179e-02
      Linear solve converged due to CONVERGED_ATOL iterations 24
    1 SNES Function norm 3.941660769344e-11
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.001 time 0.001
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

At the beginning of the output written to the terminal, we see that PyLith is reading the mesh using the `MeshIOPetsc` reader and that it found the domain to extend from 410000 to 490000 in the x direction, from 3.91e+06 to 3.99e+06 in the y direction, and from -40000 to 0 in the z direction.
The scales for nondimensionalization remain the default values for a quasistatic problem.
PyLith detects the presence of a fault based on the Lagrange multiplier for the fault in the solution field and selects appropriate preconditioning options as discussed in {ref}`sec-user-run-pylith-petsc-options`.

At the end of the output written to the terminal, we see that the solver advanced the solution one time step (static simulation).
The linear solve converged after 24 iterations and the norm of the residual met the absolute convergence tolerance (`ksp_atol`) .
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).

### Earthquake rupture parameters

We use the `pylith_eqinfo` utility to compute rupture information, such as earthquake magnitude, seismic moment, seismic potency, and average slip.
The file `eqinfoapp.cfg` holds the parameters for `pylith_eqinfo` in this example.
By default, `pylith_eqinfo` extracts information for the final time step in the output file.
In order to compute the seismic moment and moment magnitude, we need the shear modulus, which in this case is uniform over the domain.
Consquently, we specify the density and shear wave speed using a `UniformDB` in `eqinfoapp.cfg`.
The average slip, rupture area, seismic moment, and seismic potency are all given in SI units.

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
all.ruparea = [  8.140700e+08]
all.potency = [  2.623065e+09]
all.moment = [  5.901897e+19]
all.avgslip = [  3.222162e+00]
all.mommag = [  7.147328e+00]
main_fault = RuptureStats()
main_fault.timestamp = [  3.155760e+07]
main_fault.ruparea = [  5.424291e+08]
main_fault.potency = [  2.169716e+09]
main_fault.moment = [  4.881862e+19]
main_fault.avgslip = [  4.000000e+00]
main_fault.mommag = [  7.092390e+00]
east_branch = RuptureStats()
east_branch.timestamp = [  3.155760e+07]
east_branch.ruparea = [  8.993288e+07]
east_branch.potency = [  8.993288e+07]
east_branch.moment = [  2.023490e+18]
east_branch.avgslip = [  1.000000e+00]
east_branch.mommag = [  6.170734e+00]
west_branch = RuptureStats()
west_branch.timestamp = [  3.155760e+07]
west_branch.ruparea = [  1.817081e+08]
west_branch.potency = [  3.634162e+08]
west_branch.moment = [  8.176864e+18]
west_branch.avgslip = [  2.000000e+00]
west_branch.mommag = [  6.575058e+00]
```

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:crustal:strikeslip:3d:step01:solution` we use the `pylith_viz` utility to visualize the y displacement field.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step01_slip-domain.h5 warp_grid --component=y
```

:::{figure-md} fig:example:crustal:strikeslip:3d:step01:solution
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
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:41:read
 -- meshiocubit(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:185:virtual void pylith::meshio::MeshIOCubit::_read()
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Read 5253 vertices.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:189:virtual void pylith::meshio::MeshIOCubit::_read()
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Read 27596 cells.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:356:void pylith::meshio::MeshIOCubit::_readNodeSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Found 3 node sets.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:382:void pylith::meshio::MeshIOCubit::_readNodeSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading node set 'fault_main_edges' with id 30 containing 27 nodes.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:382:void pylith::meshio::MeshIOCubit::_readNodeSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading node set 'fault_west_edges' with id 31 containing 18 nodes.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:382:void pylith::meshio::MeshIOCubit::_readNodeSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading node set 'fault_east_edges' with id 32 containing 15 nodes.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:407:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Found 9 side sets.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:439:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading side set 'boundary_south' with id 10 containing 286 faces.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:439:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading side set 'boundary_east' with id 11 containing 179 faces.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:439:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading side set 'boundary_north' with id 12 containing 290 faces.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:439:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading side set 'boundary_west' with id 13 containing 172 faces.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:439:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading side set 'boundary_bottom' with id 14 containing 583 faces.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:439:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading side set 'boundary_top' with id 15 containing 930 faces.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:439:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading side set 'fault_main' with id 20 containing 197 faces.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:439:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading side set 'fault_west' with id 21 containing 72 faces.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:439:void pylith::meshio::MeshIOCubit::_readSideSets(ExodusII &)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Reading side set 'fault_east' with id 22 containing 38 faces.
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:75:void pylith::meshio::MeshIO::read(pylith::topology::Mesh *, const bool)
 -- meshiocubit(info)
 -- Component 'meshiocubit.reader': Domain bounding box:
    (413700, 493700)
    (3.917e+06, 3.977e+06)
    (-40000, 0)


-- many lines omitted --

 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/TimeDependent.py:132:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.001 time 0.
    0 SNES Function norm 2.003292919077e-02
      Linear solve converged due to CONVERGED_ATOL iterations 19
    1 SNES Function norm 1.024898741592e-10
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.001 time 0.001
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The `MeshIOCubit` reader includes diagnostic information in the journal output related to the sizes of the nodesets and material blocks.
