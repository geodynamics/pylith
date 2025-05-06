# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`, which contains numerous comments, so we only summarize the parameters here.

## Metadata, Mesh, and Output

The `pylithapp.metadata` section specifies metadata common to all simulations in the directory.
We control the verbosity of the output written to stdout using `journal.info`.
We set the parameters for importing the finite-element mesh in `pylithapp.mesh_generator`. 

## Physics

These quasi-static simulations solve the elasticity equation and include a fault, so we have a solution field with both displacement and Lagrange multiplier subfields.
%
\begin{gather}
\vec{s} = \left(\begin{array}{c} \vec{u} \quad \vec{\lambda} \end{array}\right)^T \\
\boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u}) = \vec{0}
\end{gather}
%
We use the default `TimeDependent` problem and solution field with a single displacement subfield of basis order 1.
In addition to output of the solution over the domain, we output the solution over the -y and +y boundaries.

```{code-block} cfg
---
caption: Parameters for the solution and output of the solution on the -y and +y boundaries.
---
[pylithapp.problem]
solution = pylith.problems.SolnDispLagrange

[pylithapp.problem]
solution_observers = [domain, top_boundary, bot_boundary]
solution_observers.top_boundary = pylith.meshio.OutputSolnBoundary
solution_observers.bot_boundary = pylith.meshio.OutputSolnBoundary

[pylithapp.problem.solution_observers.top_boundary]
label = boundary_ypos
label_value = 13

[pylithapp.problem.solution_observers.bot_boundary]
label = boundary_yneg
label_value = 12
```

We use the same material properties in all of the simulations in this directory, so we specify them in `pylithapp.cfg` to avoid repeating the information in the file with parameters for each simulation.
We have two materials with a contrast in the shear modulus across the fault.

```{code-block} cfg
---
caption: Material parameters common to all simulations in this directory.
---
[pylithapp.problem]
materials = [elastic_xneg, elastic_xpos]

[pylithapp.problem.materials.elastic_xneg]
label_value = 1

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Elastic properties xneg
db_auxiliary_field.values = [density, vs, vp]
db_auxiliary_field.data = [2500.0*kg/m**3, 3.00*km/s, 5.29*km/s]

auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0

derived_subfields.cauchy_strain.basis_order = 0
derived_subfields.cauchy_stress.basis_order = 0

[pylithapp.problem.materials.elastic_xpos]
label_value = 2

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Elastic properties xpos
db_auxiliary_field.values = [density, vs, vp]
db_auxiliary_field.data = [2500.0*kg/m**3, 4.24*km/s, 5.29*km/s]

auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0

derived_subfields.cauchy_strain.basis_order = 0
derived_subfields.cauchy_stress.basis_order = 0
```

Similarly, we set the general fault parameters common to all simulations in the directory.

```{code-block} cfg
---
caption: General fault parameters common to all simulations in this directory.
---
[pylithapp.problem]
interfaces = [fault]

[pylithapp.problem.interfaces.fault]
label = fault
label_value = 20

# Output `slip` and the change in fault traction on the fault.
observers.observer.data_fields = [slip, traction_change]
```

All simulations in this directory use Dirichlet boundary conditions on the -x and +x boundaries with zero displacements,
%
\begin{gather}
u_x(-50km,y) = 0,\\
u_y(-50km,y) = 0,\\
u_x(+50km,y) = 0,\\
u_y(+50km,y) = 0.
\end{gather}

```{code-block} cfg
---
caption: Dirichlet boundary conditions common to all simulations in this directory.
---
[pylithapp.problem]
bc = [bc_xneg, bc_xpos]
bc.bc_xneg = pylith.bc.DirichletTimeDependent
bc.bc_xpos = pylith.bc.DirichletTimeDependent

[pylithapp.problem.bc.bc_xpos]
label = boundary_xpos
label_value = 11
constrained_dof = [0, 1]
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC +x boundary

[pylithapp.problem.bc.bc_xneg]
label = boundary_xneg
label_value = 10
constrained_dof = [0, 1]
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC -x boundary
```
