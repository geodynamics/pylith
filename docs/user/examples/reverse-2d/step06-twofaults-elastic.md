# Step 6: Slip on Two Faults and Elastic Materials

% Metadata extracted from parameter files.
```{include} step06_twofaults_elastic-synopsis.md
```

In this example we add coseismic slip on the splay fault with an origin time of 40 years.
We specify 2 meters of reverse slip on the main fault and 1 meter of reverse slip on the splay fault.
{numref}`fig:example:reverse:2d:step06:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step06_twofaults-elastic.cfg`.

:::{figure-md} fig:example:reverse:2d:step06:diagram
<img src="figs/step06-diagram.*" alt="" scale="75%">

Boundary conditions for static coseismic slip on both the main and splay faults.
We prescribe 2 meters of reverse slip on the main fault with 1 meter of reverse slip on the splay fault.
We use roller boundary conditions on the lateral sides and bottom of the domain.
:::

## Simulation parameters

We use uniform refinement to reduce the discretization size by a factor of 2 and increase the numerical resolution.
With different origin times for slip on the two faults, we need to add time stepping.
With a linearly elastic, quasistatic model (no inertia) and fixed displacements on the boundaries, the displacement field only changes as a result of fault slip, so we can use just a few time steps to resolve the deformation.
We impose slip on the main fault in the first time step (advancing the solution to t=0), and impose slip on the splay fault in the third time step (advancing the solution to t=40 years).

```{code-block} cfg
---
caption: Parameters for uniform refinement and time stepping in Step 6.
---
[pylithapp.mesh_generator]
refiner = pylith.topology.RefineUniform

[pylithapp.problem]
initial_dt = 20.0*year
start_time = -20.0*year
end_time = 40.0*year
```

:::{important}
In 2D simulations slip is specified in terms of opening and left-lateral components.
This provides a consistent, unique sense of slip that is independent of the fault orientation.
For our geometry in this example, right lateral slip corresponds to reverse slip on both of the dipping faults.
:::

:::{important}
Both faults contain one end that is buried within the domain.
The splay fault ends where it meets the main fault.
When PyLith inserts cohesive cells into a mesh with buried edges (in this case a point), we must identify these buried edges so that PyLith properly adjusts the topology along these edges.

For properly topology of the cohesive cells, the main fault _must_ be listed first in the array of faults so that it will be created before the splay fault.
:::

We create an array of 2 faults, which are `FaultCohesiveKin` by default, and use `UniformDB` objects to specify uniform reverse slip on each fault.

```{code-block} cfg
---
caption: Parameters for prescribed earthquake rupture on the main and splay faults for Step 6.
---
[pylithapp.problem]
interfaces = [fault, splay]

[pylithapp.problem.interfaces.fault]
label = fault
label_value = 20
edge = fault_end
edge_value = 21

observers.observer.data_fields = [slip]

[pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture for main fault
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]

[pylithapp.problem.interfaces.splay]
label = splay
label_value = 22
edge = splay_end
edge_value = 23

observers.observer.data_fields = [slip]

[pylithapp.problem.interfaces.splay.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Fault rupture for splay fault
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, -1.0*m, 0.0*m]
```

## Running the simulation

```{code-block} console
---
caption: Run Step 6 simulation
---
$ pylith step06_twofaults_elastic.cfg

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
    (-100000, 100000)
    (-100000, 0)

# -- many lines omitted --

 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/TimeDependent.py:132:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 2.225574998436e-02
    Linear solve converged due to CONVERGED_ATOL iterations 41
    1 SNES Function norm 2.051112304555e-12
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

From the end of the output written to the terminal window, we see that the linear solver converged in 41 iterations and met the absolute convergence tolerance (`ksp_atol`).
As we expect for this linear problem, the nonlinear solver converged in 1 iteration.

## Visualizing the results

In {numref}`fig:example:reverse:2d:step06:solution` we use the `pylith_viz` utility to visualize the x displacement field.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filenames=output/step06_twofaults_elastic-domain.h5 warp_grid --component=x
```

:::{figure-md} fig:example:reverse:2d:step06:solution
<img src="figs/step06-solution.*" alt="Solution for Step 6 at t=40 yr. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 1000." width="600px"/>

Solution for Step 6 at t=40 yr.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is shown by the gray wireframe.
:::
