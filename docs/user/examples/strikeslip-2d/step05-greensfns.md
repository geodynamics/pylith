# Step 5: Green's Functions

% Metadata extracted from parameter files.
```{include} step05_greensfns-synopsis.md
```

## Simulation parameters

In this example we compute static Green's functions for fault slip and use then in Step 6 to invert for fault slip.
We generated the "observations" for the slip inversion in Step 4.

We impose fault slip impulses over the central portion of the strike-slip fault (-25 km $\le$ y $\le$ +25km), which is slightly larger than where we specified coseismic in Step 4. {numref}`fig:example:strikeslip:2d:step05:diagram` summarizes the boundary conditions and fault slip.
The parameters specific to this example are in `step05_greensfns.cfg`.

:::{figure-md} fig:example:strikeslip:2d:step05:diagram
<img src="figs/step05-diagram.*" alt="" scale="75%">

Boundary conditions for static Green's functions.
We set the x and y displacement to zero on the +x and -x boundaries and prescribe left-lateral slip impulses.
:::

We use the `GreensFns` problem and specify the fault on which to impose fault slip impulses.
As in Step 4, we include output at the fake GNSS stations using `OutputSolnPoints`.
In the fault interfaces section we set the fault type to `FaultCohesiveImpulses` for our fault where we want to impose fault slip impulses for the Green's functions.
We also use a spatial database to limit the section of the fault where we impose the fault slip impulses to -25 km $\le$ y $\le$ +25 km.

:::{important}
**Currently, we recommend a basis order of 1 (default) for the slip auxiliary subfield.**

The basis order for the slip auxiliary subfield controls the representation of the slip field for the impulses.
For a given impulse, a basis order of 1 will impose unit slip at a vertex with zero slip at all other vertices.
A basis order of 0 will impose unit slip over a cell with zero slip in all other cells; however, this creates jumps in the slip at the cell boundaries and reduces the accuracy of the slip inversion for a given mesh resolution.
A basis order of 2 will impose slip at vertices as well as edge degrees of freedom in the cell; regularization becomes trickier with the higher order discretization, and a basis order of 2 does not always give more accurate results for a given number of impulses.
:::

```{code-block} cfg
---
caption: Parameters for computing static Green's functions for fault slip impulses. We change the problem type and specify the fault on which we apply the slip impulses.
---
[pylithapp]
problem = pylith.problems.GreensFns

[pylithapp.greensfns]
label = fault
label_value = 20
```

```{code-block} cfg
---
caption: Parameters for the slip impulses. We change the fault type and limit the impulses to the lateral slip component.
---
[pylithapp.problem.interfaces]
fault = pylith.faults.FaultCohesiveImpulses

[pylithapp.problem.interfaces.fault]
impulse_dof = [1]

db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.iohandler.filename = slip_impulses.spatialdb

auxiliary_subfields.slip.basis_order = 1
```

## Running the simulation

```{code-block} console
---
caption: Run Step 5 simulation
---
$  pylith step05_greensfns.cfg

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

 >> /src/cig/pylith/libsrc/pylith/problems/GreensFns.cc:322:void pylith::problems::GreensFns::solve()
 -- greensfns(info)
 -- Component 'greensfns.problem': Computing Green's function 13 of 13.
  0 SNES Function norm 3.827089524877e-02
    Linear solve converged due to CONVERGED_ATOL iterations 11
  1 SNES Function norm 1.378783338243e-09
  Nonlinear solve converged due to CONVERGED_ITS iterations 1
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- greensfns(info)
 -- Finalizing problem.
WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
There are 3 unused database options. They are:
Option left: name:-ts_error_if_step_fails (no value) source: code
Option left: name:-ts_monitor (no value) source: code
Option left: name:-ts_type value: beuler source: code
```

The beginning of the output written to the terminal matches that in our previous simulations.
The second half of the output written to the terminal resembles the output from time-dependent problems, but with the time step information replaced by the impulse information.
The journal info associated with the `GreensFns` component (`journal.info.greensfns`) turns on the impulse information.
We get warnings about unused PETSc options because we do not use time stepping.

## Visualizing the results

In {numref}`fig:example:strikeslip:2d:step05:impulses` we use the `viz/plot_slip_impulses.py` Python script to visualize the slip impulses.
In {numref}`fig:example:strikeslip:2d:step05:solution` we use the `pylith_viz` utility to visualize the y displacement field.
You can move the slider or use the `p` and `n` keys to change the increment or decrement slip impulses (shown as different time stamps).

```{code-block} console
---
caption: Plot the slip impulses and visualize PyLith output using `pylith_viz`.
---
viz/plot_slip_impulses.py
pylith_viz --filename=output/step05_greensfns-domain.h5 warp_grid --component=y
```

:::{figure-md} fig:example:strikeslip:2d:step05:impulses
<img src="figs/step05_greensfns-impulses.*" alt="Slip impulses for the Green's functions in Step 5." width="400px"/>

Slip impulses for the Green's functions in Step 5.
:::

:::{figure-md} fig:example:strikeslip:2d:step05:solution
<img src="figs/step05-solution.*" alt="Solution for Step 5. The colors indicate the y displacement, and the deformation is exaggerated by a factor of 1000." width="400px"/>

Solution for Step 4.
The colors of the shaded surface indicate the y displacement, and the deformation is exaggerated by a factor of 1000.
The time value corresponds to the zero-based index of the slip impulses.
:::
