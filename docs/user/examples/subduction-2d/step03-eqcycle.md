# Step 3: Quasistatic Earthquake Cycle

% Features extracted from simulation parameter files.
```{include} step03_eqcycle-synopsis.md
```

## Simulation parameters

This simulation combines 300 years of interseismic deformation from Step 2 with the coseismic deformation from Step 1 applied at 150 years to create a simple model of an earthquake cycle.
{numref}`fig:example:subduction:2d:step03:diagram` shows the schematic of the boundary conditions.
The parameters specific to this example are in `step03_eqcycle.cfg`.

:::{figure-md} fig:example:subduction:2d:step03:diagram
<img src="figs/step03-diagram.*" alt="" width="75%">

Boundary conditions for a simple earthquake cycle with presribed coseismic slip and creep.
We combine the coseismic slip from Step 1 with the interseismic slip from Step 2.
:::

On the interface along the top of the subducting oceanic crust and the continental crust and mantle we create two earthquake ruptures.
The first rupture applies the coseismic slip from Step 1 at 150 years, while the second rupture prescribes the same steady, aseismic slip as in Step 2.
On the interface between the bottom of the subducting oceanic crust and the mantle, we prescribe the same steady, aseismic slip as that in Step 2.

```{code-block} cfg
---
caption: Prescribed slip parameters for Step 3. We only show the details for the top of the slab.
---
[pylithapp.problem]
interfaces = [fault_slabtop, fault_slabbot]

[pylithapp.problem.interfaces.fault_slabtop]
label = fault_slabtop
label_value = 21
edge = fault_slabtop_edge
edge_value = 31

observers.observer.data_fields = [slip, traction_change]

eq_ruptures = [creep, earthquake]

# Creep
eq_ruptures.creep = pylith.faults.KinSrcConstRate
eq_ruptures.creep.origin_time = 0.0*year


# Earthquake
eq_ruptures.earthquake = pylith.faults.KinSrcStep
eq_ruptures.earthquake.origin_time = 150.0*year

[pylithapp.timedependent.interfaces.fault_slabtop.eq_ruptures.earthquake]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.iohandler.filename = fault_coseismic.spatialdb
db_auxiliary_field.query_type = linear

[pylithapp.timedependent.interfaces.fault_slabtop.eq_ruptures.creep]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Fault rupture auxiliary field spatial database
db_auxiliary_field.iohandler.filename = fault_slabtop_creep.spatialdb
db_auxiliary_field.query_type = linear
```

## Running the simulation

```{code-block} console
---
caption: Run Step 3 simulation
---
$ pylith step03_eqcycle.cfg

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
    (-600000, 600000)
    (-600000, 399.651)

# -- many lines omitted --

60 TS dt 0.05 time 2.95
    0 SNES Function norm 5.748198611842e-01
      Linear solve converged due to CONVERGED_ATOL iterations 1
    1 SNES Function norm 2.523053498061e-08
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
61 TS dt 0.05 time 3.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The beginning of the output written to the terminal is identical to that from Steps 1 and 2.
At the end of the output, we see that the simulation advanced the solution 61 time steps.
Remember that the PETSc TS monitor shows the nondimensionalized time and time step values.

## Visualizing the results

In {numref}`fig:example:subduction:2d:step03:solution` we use the `pylith_viz` utility to visualize the x displacement field.
You can move the slider or use the `p` and `n` keys to change the increment or decrement time.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step03_eqcycle-domain.h5 warp_grid --component=x
```

:::{figure-md} fig:example:subduction:2d:step03:solution
<img src="figs/step03-solution.*" alt="Solution for Step 3 at t=200 yr. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 1000." width="600px"/>

Solution for Step 3 at t=200 yr.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 1000.
:::
