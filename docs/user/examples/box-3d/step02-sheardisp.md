# Step 2: Shear Displacement

% Meatadata extracted from parameter files
```{include} step02_sheardisp-synopsis.md
```

## Simulation parameters

This example corresponds to shear deformation due to Dirichlet (displacement) boundary conditions.
We apply Dirichlet (displacement) boundary conditions for the y displacement on the +x (`boundary_xpos`) and -x (`boundary_xneg`) boundaries and for the x displacement on the +y (`boundary_ypos`) and -y (`boundary_yneg`) boundaries.
{numref}`fig:example:box:3d:step01:diagram` shows the boundary conditions on the domain.
The parameters specific to this example are in `step02_sheardisp.cfg`.

:::{figure-md} fig:example:box:3d:step02:diagram
<img src="figs/step02-diagram.*" alt="" scale="75%">

Boundary conditions for shear deformation.
We constrain the y displacement on the +x and -x boundaries and the x displacement on the +y and -y boundaries.
:::

We create an array of 5 Dirichlet boundary conditions.
On the -x, +x, -y, and +y boundaries we impose shear displacement.

```{code-block} cfg
---
caption: Boundary conditions for Step 2. We only show the details for the -x boundary.
---
bc = [bc_xneg, bc_xpos, bc_yneg, bc_ypos, bc_zneg]
bc.bc_xneg = pylith.bc.DirichletTimeDependent
bc.bc_xpos = pylith.bc.DirichletTimeDependent
bc.bc_yneg = pylith.bc.DirichletTimeDependent
bc.bc_ypos = pylith.bc.DirichletTimeDependent
bc.bc_zneg = pylith.bc.DirichletTimeDependent

[pylithapp.problem.bc.bc_xneg]
label = boundary_xneg
label_value = 10
constrained_dof = [1]

db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Dirichlet BC -x boundary
db_auxiliary_field.iohandler.filename = sheardisp_bc_xneg.spatialdb
db_auxiliary_field.query_type = linear
```

## Running the simulation

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02_sheardisp.cfg

# The output should look something like the following.
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/meshio/MeshIOObj.py:44:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:94:void pylith::meshio::MeshIO::read(topology::Mesh *)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-6000, 6000)
    (-6000, 6000)
    (-9000, 0)

# -- many lines omitted --

 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/TimeDependent.py:139:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 1.811215061775e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 1
    1 SNES Function norm 2.330640615892e-17 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

The output written to the terminal is nearly identical to what we saw for Step 1.

## Visualizing the results

In {numref}`fig:example:box:3d:step02:solution` we use the `pylith_viz` utility to visualize the x displacement field.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filenames=output/step02_sheardisp-domain.h5 warp_grid --component=x
```

:::{figure-md} fig:example:box:3d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 1000." width="600px"/>

Solution for Step 2.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is show by the gray wireframe.
:::
