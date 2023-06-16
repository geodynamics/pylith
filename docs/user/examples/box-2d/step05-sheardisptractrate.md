# Step 5: Time-Dependent Shear Displacement and Tractions

% Meatadata extracted from parameter files
```{include} step05_sheardisptractrate-synopsis.md
```

## Simulation parameters

In this example we build on Step 3 and make the Dirichlet (displacement) and Neumann (traction) boundary conditions a bit more complicated by adding variation in time.
The simulation has a duration of 5 years with a time step of 1 year.
The time-dependent boundary conditions use the same initial amplitude values for the first time step before adding in a constant rate increase at a time of 1 year.
{numref}`fig:example:box:2d:step05:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step05_sheardisptractrate.cfg`.

:::{figure-md} fig:example:box:2d:step05:diagram
<img src="figs/step05-diagram.*" alt="" scale="75%">

Boundary conditions for shear deformation.
We constrain the x and y displacements on the +x and -x boundaries.
We apply tangential (shear) tractions on the +y and -y boundaries.
At a time of 1 year we increase the amplitude at a constrant rate $b$ ($H(t)$ corresponds to the heavyside step function).
:::

This is a time-dependent problem, so we must specify the start and end times of the simulation along with the initial time step.
With an initial time step of 1 year, we start the simulation at -1 year so that the first solve will advance the simulation to a time of 0.
We also specify a relaxation time on the order of the time scale of the simulation to allow for reasonable nondimensionalization of time.

```{code-block} cfg
---
caption: Time stepping parameters for Step 5.
---
[pylithapp.problem]
start_time = -1.0*year
end_time = 5.0*year
initial_dt = 1.0*year

[pylithapp.problem.normalizer]
relaxation_time = 10.0*year
```

For the time-dependent Dirichlet and Neumann boundary conditions, we specify both the initial displacement and a constant rate; the constant rate begins at t=1 year.

```{code-block} cfg
---
caption: Time-dependent boundary conditions for Step 5. We show the details for the -x and -y boundaries.
---
[pylithapp.problem.bc.bc_xneg]
# Degrees of freedom (dof) 0 and 1 correspond to the x and y displacements. 
constrained_dof = [0, 1]
label = boundary_xneg
use_initial = True
use_rate = True
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Dirichlet BC -x boundary
db_auxiliary_field.iohandler.filename = sheardisprate_bc_xneg.spatialdb
db_auxiliary_field.query_type = linear

[pylithapp.problem.bc.bc_yneg]
label = boundary_yneg
use_initial = True
use_rate = True
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Neumann BC -y boundary
db_auxiliary_field.values = [initial_amplitude_tangential, initial_amplitude_normal, rate_start_time, rate_amplitude_tangential, rate_amplitude_normal]
db_auxiliary_field.data = [-4.5*MPa, 0.0*MPa, 1.0*year, -1.125*MPa/year, 0.0]
```

## Running the simulation

```{code-block} console
---
caption: Run Step 5 simulation
---
$ pylith step05_sheardisptractrate.cfg

# The output should look something like the following.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshioascii(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(topology::Mesh *)
 -- meshioascii(info)
 -- Component 'reader': Domain bounding box:
    (-6000, 6000)
    (-16000, -0)

# -- many lines omitted --

5 TS dt 0.1 time 0.4
    0 SNES Function norm 1.467261021331e-03 
    Linear solve converged due to CONVERGED_ATOL iterations 1
    1 SNES Function norm 7.893110957891e-19 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
6 TS dt 0.1 time 0.5
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The output written to the terminal now contains multiple time steps.
The PETSc TS (time stepping) monitor shows the time step number followed by the time step and time in nondimensional units.

## Visualizing the results

In {numref}`fig:example:box:2d:step03:solution` we use ParaView to visualize the x displacement field using the `viz/plot_dispwarp.py` Python script.
As in Step 2 we override the default name of the simulation file with the name of the current simulation before running the `viz/plot_dispwarp.py` Python script.

```{code-block} python
---
caption: Set the simulation in the ParaView Python Shell.
---
>>> SIM = "step05_sheardisptractrate"
```

One you run the `viz/plot_dispwarp.py` Python script, you can click on the "play" button corresponding to the right triangle in the toolbar to view the time-dependent deformation.

:::{figure-md} fig:example:box:2d:step05:solution
<img src="figs/step05-solution.*" alt="Solution for Step 5. The colors indicate the magnitude of the displacement, and the deformation is exaggerated by a factor of 1000." width="100%"/>

Solution for Step 5 at a time of 4.0 years.
The colors of the shaded surface indicate the magnitude of the x displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is show by the gray wireframe.
:::
