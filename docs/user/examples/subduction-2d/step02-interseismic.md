# Step 2: Quasistatic Interseismic Deformation

In this example we simulate the interseismic deformation associated with the oceanic crust subducting beneath the continental crust and into the mantle.
We prescribe steady aseismic slip of 8 cm/yr along the interfaces between the oceanic crust and mantle with the interface between the oceanic crust and continental crust locked as shown in {numref}`fig:example:subduction:2d:step02:diagram`.
We adjust the Dirichlet (displacement) boundary conditions on the lateral edges and bottom of the domain by pinning only the portions of the boundaries that are mantle and continental crust and not oceanic crust.

:::{figure-md} fig:example:subduction:2d:step02:diagram
<img src="figs/step02-diagram.*" alt="" width="100%">

Boundary conditions for quasistatic simulation for interseismic deformation.
We prescribe constant creep on the top and bottom of the subduction slab, except for the portion of the subduction interface where we imposed coseismic slip in Step 1.
We lock (zero creep) that part of the interface.
:::

% Features extracted from simulation parameter files.
```{include} step02_interseismic-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step02_interseismic.cfg`.
These include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters defining the start time and end time for the quasistatic simulation.
* `pylithapp.problem.fault` Parameters for prescribed slip on the fault.
* `pylithapp.problem.bc` Parameters for velocity boundary conditions.

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02_interseismic.cfg

# The output should look something like the following.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(topology::Mesh *)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-600000, 600000)
    (-600000, 399.651)

# -- many lines omitted --

30 TS dt 0.05 time 1.45
    0 SNES Function norm 5.748198604376e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 178
    1 SNES Function norm 1.124343852602e-11 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
31 TS dt 0.05 time 1.5
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The beginning of the output written to the terminal is identical to that from Step 1.
At the end of the output, we see that the simulation advanced the solution 31 time steps.
Remember that the PETSc TS monitor shows the nondimensionalized time and time step values.

## Visualizing the results

In {numref}`fig:example:subduction:2d:step02:solution` we use ParaView to visualize the x displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/subduction-2d` directory.
Next, we override the default name of the simulation file with the name of the current simulation.

```{code-block} python
---
caption: Set the simulation in the ParaView Python Shell.
---
>>> SIM = "step02_interseismic"
```

Finally, we run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.

:::{figure-md} fig:example:subduction:2d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2 at t=100 yr. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 2 at t=100 yr.
The colors of the shaded surface indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000.
:::
