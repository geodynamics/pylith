# Step 2: No faults, with flexure

```{include} step02-no-faults-with-flexure-synopsis.md
```

## Simulation parameters

This example uses poroelasticity to model the infiltration of seawater through a slab of oceanic lithosphere. The permeability field is depth dependent, decreasing with depth but does not vary laterally. The lithosphere is now subject to deformation, over 300 kyr the slab bends to simulate extensional stresses in the outer-rise of a subduction zone. A fluid pressure is applied to the top boundary that is equivalent to the pressure exerted on the seafloor by the water column. This simulates what the hydration state of the oceanic lithosphere as it is about to enter a convergent margin. 


{numref}`fig:example:poroelastic:outerrise:2d:step02:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step02-no-faults-with-flexure.cfg`.

:::{figure-md} fig:example:outerrise:poroelastic:2d:step02:diagram
<img src="figs/step02-diagram.*" alt="" scale="100%">

Boundary and initial conditions for step02.
We fix the left boundary, but we now apply a spatially varying velocity condition on the top boundary using a `SimpleDB` file, while leaving the right and bottom boundaries unconstrained. 
We impose a fluid pressure  on the +y boundary equal to the weight of the water column to generate fluid flow.
:::

```{code-block} cfg
---
caption: Modified top boundary condition for Step 2.
---
[pylithapp.problem.bc.bndry_top]
use_initial = False
use_rate = True

db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Dirichlet BC +y boundary
db_auxiliary_field.iohandler.filename = top_velocity_boundary.spatialdb
```

```{code-block} cfg
---
caption: Time stepping parameters for Step 2.
---
[pylithapp.timedependent]
start_time = -6e3*year
initial_dt = 6e3*year
end_time = 300e3*year
```

```{code-block} cfg
---
caption: Initial condition parameters for Step 2. We use a `SimpleGridDB` file that does not contain enhanced permeability due to outer rise faults. 
---
[pylithapp.problem]

[pylithapp.problem.materials.slab]
db_auxiliary_field.filename = no_faultzone_permeability.spatialgrid
```

## Running the simulation

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02-no-faults-with-flexure.cfg

 software/pylith-debug/lib/python3.11/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 software/pylith-debug/lib/python3.11/site-packages/pylith/meshio/MeshIOObj.py:38:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /pylith-main/libsrc/pylith/meshio/MeshIO.cc:85:void pylith::meshio::MeshIO::read(pylith::topology::Mesh *, const bool)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (0, 150000)
    (-30000, 0)

# -- many lines omitted --

>> /software/pylith-debug/lib/python3.11/site-packages/pylith/problems/TimeDependent.py:132:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 6000. time -6000.
    0 SNES Function norm 3.746797722103e+02 
    Linear solve converged due to CONVERGED_RTOL iterations 25
    1 SNES Function norm 2.600147960846e-10 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 6000. time 0.
    0 SNES Function norm 1.697747712084e+02 
    Linear solve converged due to CONVERGED_RTOL iterations 23
    1 SNES Function norm 5.641202104179e-10 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1

# -- many lines omitted --

50 TS dt 6000. time 294000.
    0 SNES Function norm 1.697753481522e+02 
    Linear solve converged due to CONVERGED_RTOL iterations 14
    1 SNES Function norm 8.568690135687e-10 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
51 TS dt 6000. time 300000.
 >> /software/pylith-debug/lib/python3.11/site-packages/pylith/problems/Problem.py:199:finalize
 -- timedependent(info)
 -- Finalizing problem.
```


:::{figure-md} fig:example:poroelastic:outerrise:2d:gmsh:mesh
<img src="figs/step02-porosity.*" alt="Porosity field at the end of model run time for step02." width="100%"/>

Porosity field at the end of model run time for step02.
:::