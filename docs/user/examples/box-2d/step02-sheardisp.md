# Step 2: Shear Displacement

This example corresponds to shear deformation due to Dirichlet (displacement) boundary conditions.
We apply Dirichlet (displacement) boundary conditions for the y displacement on the +x (`boundary_xpos`) and -x (`boundary_xneg`) boundaries and for the x displacement on the +y (`boundary_ypos`) and -y (`boundary_yneg`) boundaries.
{numref}`fig:example:box:2d:step01:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:box:2d:step02:diagram
<img src="figs/step02-diagram.*" alt="" scale="75%">

Boundary conditions for shear deformation.
We constrain the y displacement on the +x and -x boundaries and the x displacement on the +y and -y boundaries.
:::

% Meatadata extracted from parameter files
```{include} step02_sheardisp-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step02_sheardisp.cfg`.
These include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem.solution` Specify the basis order for the solution fields, in this case the `displacement` field.
* `pylithapp.problem.bc` Parameters for the boundary conditions. The displacement field varies  along the boundary, so we use a `SimpleDB` spatial database and the `linear` query type.

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02_sheardisp.cfg

# The output should look something like the following.
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshioascii(info)
 -- Reading finite-element mesh
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(topology::Mesh *)
 -- meshioascii(info)
 -- Component 'reader': Domain bounding box:
    (-6000, 6000)
    (-16000, -0)

# -- many lines omitted --

 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 2.239977678460e-03 
    Linear solve converged due to CONVERGED_ATOL iterations 1
    1 SNES Function norm 1.964321818484e-18 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The output written to the terminal is nearly identical to what we saw for Step 1.
We omit the middle portion of the output which shows that the domain, the scales for nondimensionalization, and PETSc options all remain the same.

## Visualizing the results

In {numref}`fig:example:box:2d:step02:solution` we use ParaView to visualize the x displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/box-2d` directory.

```{code-block} console
---
caption: Open ParaView using the command line.
---
$ PATH_TO_PARAVIEW/paraview

# For macOS, it will be something like
$ /Applications/ParaView-5.9.1.app/Contents/MacOS/paraview
```

Next, we override the default name of the simulation file with the name of the current simulation.

```{code-block} python
---
caption: Set the simulation in the ParaView Python Shell.
---
>>> SIM = "step02_sheardisp"
```

Finally, we run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.

:::{figure-md} fig:example:box:2d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 2.
The colors of the shaded surface indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is show by the gray wireframe.
:::
