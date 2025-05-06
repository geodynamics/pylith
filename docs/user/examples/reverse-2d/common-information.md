# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`, which contains numerous comments, so we only summarize the parameters here.

We output the solution over the domain and the ground surface (+y boundary).

```{code-block} cfg
---
caption: Parameters for output of the solution over the domain and ground surface (+y boundary).
---
[pylithapp.problem]
solution_observers = [domain, boundary]
solution_observers.boundary = pylith.meshio.OutputSolnBoundary

[pylithapp.problem.solution_observers.boundary]
label = boundary_ypos
label_value = 13
```

These static and quasi-static simulations solve the elasticity equation.
We use the same material properties for several simulations in this directory, so we specify them in `pylithapp.cfg` to avoid repeating the information in the file with parameters for each simulation.
We use a `SimpleDB` spatial database so that we can simply use a different spatial database file when we change the bulk rheology.

```{code-block} cfg
---
caption: Material parameters for isotropic, linear elasticity. We only show the details for the slab material.
---
[pylithapp.problem]
materials = [slab, crust, wedge]

[pylithapp.problem.materials.slab]
label_value = 1

db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Elastic properties for slab
db_auxiliary_field.iohandler.filename = mat_elastic.spatialdb

auxiliary_subfields.density.basis_order = 0
auxiliary_subfields.gravitational_acceleration.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0
```

Similarly, for all of the simulations in this directory we use Dirichlet (displacement) boundary conditions on the +x, -x, and -y boundaries that constrain the displacement component perpendicular to the fault.

```{code-block} cfg
---
caption: Dirichlet boundary condition parameters common to all simulations in this directory. We only show the details for the +x boundary.
---
[pylithapp.problem]
bc = [bc_xneg, bc_xpos, bc_yneg]
bc.bc_xneg = pylith.bc.DirichletTimeDependent
bc.bc_xpos = pylith.bc.DirichletTimeDependent
bc.bc_yneg = pylith.bc.DirichletTimeDependent

[pylithapp.problem.bc.bc_xpos]
label = boundary_xpos
label_value = 11
constrained_dof = [0]
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC +x edge

auxiliary_subfields.initial_amplitude.basis_order = 0
```
