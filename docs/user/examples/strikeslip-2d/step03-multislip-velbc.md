# Step 3: Multiple Earthquake Ruptures and Velocity Boundary Conditions

% Metadata extracted from parameter files.
```{include} step03_multislip_velbc-synopsis.md
```

## Simulation parameters

This example involves a quasistatic simulation that solves for the deformation from velocity boundary conditions and multiple earthquake ruptures on the fault.
The velocity boundary conditions match those in Step 2.
We prescribe the first earthquake rupture to occur at 100 years with 1 meter of right-lateral slip and the second earthquake rupture to occur at 200 years with 3 meters of right-lateral slip.
{numref}`fig:example:strikeslip:2d:step03:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step03_multislip_velbc.cfg`.

:::{figure-md} fig:example:strikeslip:2d:step03:diagram
<img src="figs/step03-diagram.*" alt="" scale="75%">

Boundary conditions for quasi-static simulation with velocity boundary conditions and coseismic slip.
We set the x displacement to zero on the +x and -x boundaries.
We set the y velocity to -1 cm/yr on the +x boundary and +1 cm/yr on the -x boundary.
We prescribe 1 meter of right-lateral slip to occur at 100 years and 3 meters of right-lateral slip to occur at 200 years.
:::

```{code-block} cfg
---
caption: Prescribed slip rupture parameters for Step 3 with two earthquake ruptures on the strike-slip fault.
---
[pylithapp.problem.interfaces.fault]
eq_ruptures = [one, two]

[pylithapp.problem.interfaces.fault.eq_ruptures.one]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture one
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [100.0*year, -1.0*m, 0.0*m]

[pylithapp.problem.interfaces.fault.eq_ruptures.two]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture two
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [200.0*year, -3.0*m, 0.0*m]
```

## Running the simulation

```{code-block} console
---
caption: Run Step 3 simulation
---
$ pylith step03_multislip_velbc.cfg

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
    (-50000, 50000)
    (-75000, 75000)
# -- many lines omitted --

25 TS dt 0.1 time 2.4
    0 SNES Function norm 3.175726348635e+00
      Linear solve converged due to CONVERGED_ATOL iterations 0
    1 SNES Function norm 1.889182885363e-09
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
26 TS dt 0.1 time 2.5
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The beginning of the output written to the terminal is identical to that from Step 1.
At the end of the output, we see that the simulation advanced the solution 26 time steps.
Remember that the PETSc TS monitor shows the nondimensionalized time and time step values.

## Visualizing the results

In {numref}`fig:example:strikeslip:2d:step03:solution` we use the `pylith_viz` utility to visualize the y displacement field.
You can move the slider or use the `p` and `n` keys to change the increment or decrement time.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step03_multislip_velbc-domain.h5 warp_grid --component=y
```

:::{figure-md} fig:example:strikeslip:2d:step03:solution
<img src="figs/step03-solution.*" alt="Solution for Step 3 at t=190 yr. The colors indicate the y displacement, and the deformation is exaggerated by a factor of 1000." width="400px"/>

Solution for Step 3 at t=190 yr.
The colors of the shaded surface indicate the y displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is shown by the gray wireframe.
:::
