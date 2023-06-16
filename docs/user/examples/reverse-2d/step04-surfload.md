# Step 4: Surface Tractions

% Metadata extracted from parameter files.
```{include} step04_surfload-synopsis.md
```

## Simulation parameters

This example focuses on loading via surface tractions on the +y boundary.
We apply tractions normal to the boundary with a trapezoidal distribution as shown in {numref}`fig:example:reverse:2d:step04:diagram`
We use the same roller boundary conditions that we used in Steps 1-3.
The parameters specific to this example are in `step04_surfload.cfg`.

:::{figure-md} fig:example:reverse:2d:step04:diagram
<img src="figs/step04-diagram.*" alt="" scale="75%">

We add a Neumann (traction) boundary condition on the +y boundary with roller boundary conditions on the lateral sides and bottom of the domain.
:::

```{code-block} cfg
---
caption: Surface load parameters for Step 4.
---
[pylithapp.problem]
bc = [bc_xneg, bc_xpos, bc_yneg, bc_ypos]
bc.bc_ypos = pylith.bc.NeumannTimeDependent

[pylithapp.problem.bc.bc_ypos]
label = boundary_ypos
label_value = 13

db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Neumann BC +y edge
db_auxiliary_field.iohandler.filename = traction_surfload.spatialdb

db_auxiliary_field.query_type = linear

auxiliary_subfields.initial_amplitude.basis_order = 1
```

## Running the simulation

```{code-block} console
---
caption: Run Step 4 simulation
---
$ pylith step04_surfload.cfg

# The output should look something like the following.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(pylith::topology::Mesh*)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-100000, 100000)
    (-100000, 0)

# -- many lines omitted --

 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 1.213351093160e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 1
    1 SNES Function norm 1.038106792811e-15 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
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
