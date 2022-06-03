# Step 3: Quasistatic Earthquake Cycle

This simulation combines 300 years of interseismic deformation from Step 2 with the coseismic deformation from Step 1 applied at 150 years to create a simple model of an earthquake cycle.
{numref}`fig:example:subduction:2d:step03:diagram` shows the schematic of the boundary conditions.

:::{figure-md} fig:example:subduction:2d:step03:diagram
<img src="figs/step03-diagram.*" alt="" width="100%">

Boundary conditions for a simple earthquake cycle with presribed coseismic slip and creep.
We combine the coseismic slip from Step 1 with the interseismic slip from Step 2.
:::

% Features extracted from simulation parameter files.
```{include} step03_eqcycle-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step03_eqcycle.cfg`.
These include:

* `pylithapp.metadata` Metadata for this simulation. Even when the author and version are the same for all simulations in a directory, we prefer to keep that metadata in each simulation file as a reminder to keep it up-to-date for each simulation.
* `pylithapp` Parameters defining where to write the output.
* `pylithapp.problem` Parameters defining the start time and end time for the quasistatic simulation.
* `pylithapp.problem.fault` Parameters for prescribed slip on the fault.
* `pylithapp.problem.bc` Parameters for velocity boundary conditions.

On the interface along the top of the subducting oceanic crust and the continental crust and mantle we create two earthquake ruptures.
The first rupture applies the coseismic slip from Step 1 at 150 years, while the second rupture prescribes the same steady, aseismic slip as in Step 2.
On the interface between the bottom of the subducting oceanic crust and the mantle, we prescribe the same steady, aseismic slip as that in Step 2.

```{code-block} console
---
caption: Run Step 3 simulation
---
$ pylith step03_eqcycle.cfg

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

61 TS dt 0.05 time 3.
    0 SNES Function norm 5.748198604376e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 178
    1 SNES Function norm 1.127019123456e-11 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
62 TS dt 0.05 time 3.05
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The beginning of the output written to the terminal is identical to that from Steps 1 and 2.
At the end of the output, we see that the simulation advanced the solution 62 time steps.
Remember that the PETSc TS monitor shows the nondimensionalized time and time step values.

## Visualizing the results

In {numref}`fig:example:subduction:2d:step03:solution` we use ParaView to visualize the x displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/subduction-2d` directory.
Next, we override the default name of the simulation file with the name of the current simulation.

```{code-block} python
---
caption: Set the simulation in the ParaView Python Shell.
---
>>> SIM = "step03_eqcycle"
```

Finally, we run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.

:::{figure-md} fig:example:subduction:2d:step03:solution
<img src="figs/step03-solution.*" alt="Solution for Step 3 at t=200 yr. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 3 at t=200 yr.
The colors of the shaded surface indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000.
:::
