# Step 1: Static Coseismic Slip

This example involves a static simulation that solves for the deformation from prescribed coseismic slip on the fault.
We specify 2 meters of right-lateral slip.
{numref}`fig:example:strikeslip:2d:step01:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:strikeslip:2d:step01:diagram
<img src="figs/step01-diagram.*" alt="" scale="75%">

Boundary conditions for static coseismic slip.
We set the x and y displacement to zero on the +x and -x boundaries and prescribe 2 meters of right-lateral slip.
:::

## Features

```{include} step01-slip-features.md
```

## Simulation parameters

The parameters specific to this example are in `step01_slip.cfg`.
These include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem.fault` Parameters for prescribed slip on the fault.

```{code-block} console
---
caption: Run Step 1 simulation
---
$ pylith step01_slip.cfg

# The output should look something like the following.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(topology::Mesh *)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-50000, 50000)
    (-75000, 75000)

# -- many lines omitted --

 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 4.895713226482e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 36
        Line search: Using full step: fnorm 4.895713226482e-02 gnorm 1.198850635632e-12
    1 SNES Function norm 1.198850635632e-12 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:196:finalize
 -- timedependent(info)
 -- Finalizing problem.
 ```

At the beginning of the output written to the terminal, we see that PyLith is reading the mesh using the `MeshIOPetsc` reader and that it found the domain to extend from -50,000 m to +50,000 m in the x direction and from -75,000 m to +75,000 m in the y direction.
At the end of the output written to the termial, we see that the solver advanced the solution one time step (static simulation).
The linear solve converged after 36 iterations and the norm of the residual met the absolute convergence tolerance (`ksp_atol`) .
The nonlinear solve converged in 1 iteration, which we expect because this is a linear problem, and the residual met the absolute convergence tolerance (`snes_atol`).

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:strikeslip:2d:step01:solution` we use ParaView to visualize the y displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/strikeslip-2d` directory.

```{code-block} console
---
caption: Open ParaView using the command line.
---
$ PATH_TO_PARAVIEW/paraview

# For macOS, it will be something like
$ /Applications/ParaView-5.9.1.app/Contents/MacOS/paraview
```

Next we run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.
For Step 1 we do not need to change any of the default values.

:::{figure-md} fig:example:strikeslip:2d:step01:solution
<img src="figs/step01-solution.*" alt="Solution for Step 1. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 1.
The colors of the shaded surface indicate the magnitude of the y displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is show by the gray wireframe.
The contrast in material properties across the faults causes the asymmetry in the y displacement field.
:::

## Step 1 with Cubit Mesh

Using the Cubit mesh rather than the Gmsh mesh involves two changes:

1. Use the `MeshIOCubit` reader instead of the `MeshIOPetsc` reader and change the filename of the mesh file.
2. Set the `label_value` to 1 for boundary conditions and faults.\
   We must override the nondefault `label_value` settings in `pylithapp.cfg` that were appropriate for our Gmsh reader but are incorrect for the Cubit reader.

The file `step01_slip_cubit.cfg` provides these changes and updates the names for output.

```{code-block} console
---
caption: Run Step 1 simulation with the Cubit mesh
---
$ pylith step01_slip_cubit.cfg

# The output should look something like the following.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiocubit(info)
 -- Reading finite-element mesh
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:157:void pylith::meshio::MeshIOCubit::_readVertices(pylith::meshio::ExodusII &, pylith::scalar_array *, int *, int *) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 682 vertices.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:217:void pylith::meshio::MeshIOCubit::_readCells(pylith::meshio::ExodusII &, pylith::int_array *, pylith::int_array *, int *, int *) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 1276 cells in 2 blocks.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:279:void pylith::meshio::MeshIOCubit::_readGroups(pylith::meshio::ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Found 5 node sets.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:305:void pylith::meshio::MeshIOCubit::_readGroups(pylith::meshio::ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'fault' with id 10 containing 39 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:305:void pylith::meshio::MeshIOCubit::_readGroups(pylith::meshio::ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'boundary_xpos' with id 21 containing 24 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:305:void pylith::meshio::MeshIOCubit::_readGroups(pylith::meshio::ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'boundary_xneg' with id 22 containing 24 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:305:void pylith::meshio::MeshIOCubit::_readGroups(pylith::meshio::ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'boundary_ypos' with id 23 containing 21 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:305:void pylith::meshio::MeshIOCubit::_readGroups(pylith::meshio::ExodusII &)
 -- meshiocubit(info)
 -- Component 'reader': Reading node set 'boundary_yneg' with id 24 containing 21 nodes.
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(topology::Mesh *)
 -- meshiocubit(info)
 -- Component 'reader': Domain bounding box:
    (-50000, 50000)
    (-75000, 75000)

# -- many lines omitted --

 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 4.834519229177e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 37
        Line search: Using full step: fnorm 4.834519229177e-02 gnorm 1.046009457852e-12
    1 SNES Function norm 1.046009457852e-12 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:196:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The `MeshIOCubit` reader includes diagnostic information in the journal output related to the sizes of the nodesets and material blocks.
