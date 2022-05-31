# Step 4: Surface Tractions

This example focuses on loading via surface tractions on the +y bondary.
We apply tractions normal to the boundary with a trapezoidal distribution as shown in {numref}`fig:example:reverse:2d:step04:diagram`
We use the same roller boundary conditions that we used in Steps 1-3.

:::{figure-md} fig:example:reverse:2d:step04:diagram
<img src="figs/step04-diagram.*" alt="" scale="75%">

We add a Neumann (traction) boundary condition on the +y boundary with roller boundary conditions on the lateral sides and bottom of the domain.
:::

% Metadata extracted from parameter files.
```{include} step04_surfload-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step04_surfload.cfg`.
These include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters for specifying the basis order.
* `pylithapp.problem.bc` Add a Neumann (traction) boundary condition on the +y boundary.

```{code-block} console
---
caption: Run Step 4 simulation
---
$ pylith step04_surfload.cfg

# The output should look something like the following.
 >> /software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /home/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(pylith::topology::Mesh*)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-100000, 100000)
    (-100000, 0)

# -- many lines omitted --

 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 1.213351093160e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 1
    1 SNES Function norm 1.038106792811e-15 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

As expected from the use of the LU preconditioner and linear problem, both the linear and nonlinear solvers converged in 1 iterations.

## Visualizing the results

The `output` directory contains the simulation output.
Each "observer" writes its own set of files, so the solution over the domain is in one set of files, the boundary condition information is in another set of files, and the material information is in yet another set of files.
The HDF5 (`.h5`) files contain the mesh geometry and topology information along with the solution fields.
The Xdmf (`.xmf`) files contain metadata that allow visualization tools like ParaView to know where to find the information in the HDF5 files.
To visualize the data using ParaView or Visit, load the Xdmf files.

In {numref}`fig:example:reverse:2d:step04:solution` we use ParaView to visualize the displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/reverse-2d` directory.
Before running the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`, we set the simulation name in the ParaView Python Shell.

```{code-block} python
---
caption: Set the simulation and exaggeration in the ParaView Python Shell.
---
>>> SIM = "step04_surfload"
>>> WARP_SCALE = 500
```

:::{figure-md} fig:example:reverse:2d:step04:solution
<img src="figs/step04-solution.*" alt="Solution for Step 4. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 500." width="100%"/>

Solution for Step 4.
The colors of the shaded surface indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 500.
The undeformed configuration is show by the gray wireframe.
:::
