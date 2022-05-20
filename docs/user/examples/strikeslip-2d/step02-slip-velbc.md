# Step 2: Single Earthquake Rupture and Velocity Boundary Conditions

This example involves a quasistatic simulation that solves for the deformation from velocity boundary conditions and prescribed coseismic slip on the fault.
We let strain accumulate due to the motion of the boundaries and then release the strain by prescribing 2 meters of right-lateral slip at t=100 years.
{numref}`fig:example:strikeslip:2d:step02:diagram` shows the boundary conditions on the domain.

:::{figure-md} fig:example:strikeslip:2d:step02:diagram
<img src="figs/step02-diagram.*" alt="" scale="75%">

Boundary conditions for quasistatic simulation with velocity boundary conditions and coseismic slip.
We set the x displacement to zero on the +x and -x boundaries.
We set the y velocity to -1 cm/yr on the +x boundary and +1 cm/yr on the -x boundary.
We prescribe 2 meters of right-lateral slip to occur at 100 years to release the accumulated strain energy.
:::

% Metadata extracted from parameter files.
```{include} step02_slip_velbc-synopsis.md
```

## Simulation parameters

The parameters specific to this example are in `step02_slip_velbc.cfg`.
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
$ pylith step02_slip_velbc.cfg

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

24 TS dt 0.05 time 1.15
    0 SNES Function norm 5.390420823600e-04 
    Linear solve converged due to CONVERGED_ATOL iterations 32
    1 SNES Function norm 1.481819020131e-12 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
25 TS dt 0.05 time 1.2
 >> /Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The beginning of the output written to the terminal is identical to that from Step 1.
At the end of the output, we see that the simulation advanced the solution 25 time steps.
Remember that the PETSc TS monitor shows the nondimensionalized time and time step values.

## Visualizing the results

In {numref}`fig:example:strikeslip:2d:step02:solution` we use ParaView to visualize the x displacement field using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/strikeslip-2d` directory.

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
>>> SIM = "step02_slip_velbc"
```

Finally, we run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.

:::{tip}
You can use the "play" button to animate the solution in time.
:::

:::{figure-md} fig:example:strikeslip:2d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2 at t=100 yr. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 2 at t=100 yr.
The colors of the shaded surface indicate the magnitude of the y displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is show by the gray wireframe.
The coseismic fault slip at 100 years releases all of the accumulated strain energy.
:::
