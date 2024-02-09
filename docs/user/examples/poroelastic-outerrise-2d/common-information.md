# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`.
The `pylithapp.cfg` file contains numerous comments, so we only summarize the parameters here.

## Metadata, Mesh, and Output

The `pylithapp.metadata` section specifies metadata common to all simulations in the directory.
We control the verbosity of the output written to stdout using `journal.info`.
We set the parameters for importing the finite-element mesh in `pylithapp.mesh_generator`. 

## Physics

These quasi-static simulations solve the poroelasticity equation with state-variables, so we have a solution field with displacement, pressure, volumetric strain, velocity, the time derivative of pressure, and the volumetric strain rate subfields.
%
\begin{gather}
\vec{s} = \left(\vec{u} \quad p \quad \epsilon_v \quad \vec{v} \quad \dot{p} \quad \dot{\epsilon_v}\right)^T, \\
\nabla \cdot \boldsymbol{\sigma}(\vec{u},p) = \vec{0}, \\
\frac{\partial \zeta(\vec{u},p)}{\partial t} + \nabla \cdot \vec{q}(p) = 0, \\
\nabla \cdot \vec{u} - \epsilon_{v} = 0.
\frac{\partial \phi (\vec{x}, t)}{\partial t} = (\alpha (\vec{x}) - \phi (\vec{x}, t)) \left(\dot{\epsilon_v} + \frac{1 - \alpha (\vec{x})}{K_d(\vec{x})} \dot{p} (\vec{x}, t) \right)
\end{gather}

We specify a basis order of 2 for the displacement, the velocity, and the isotropic permeability subfields, and 1 for each of the remaining subfields.

```{code-block} cfg
---
caption: Solution and nondimensionalization parameters for poroelastic-outerrise-2d examples.
---
[pylithapp.problem]
solution = pylith.problems.SolnDispPresTracStrainVelPdotTdot
defaults.quadrature_order = 2

[pylithapp.problem.solution.subfields]
displacement.basis_order = 2
pressure.basis_order = 1
trace_strain.basis_order = 1

velocity.basis_order = 2
pressure_t.basis_order = 1
trace_strain_t.basis_order = 1

[pylithapp.problem.materials.slab.bulk_rheology]
auxiliary_subfields.isotropic_permeability.basis_order = 2
```

We use two different SimpleGridDB files to initialize the material properties for the three steps. So set all other relevant material parameters in the `pylithapp.cfg` file to avoid repeat lines. 

We only have a single material; with spatially varying material properties so we need to use a `SimpleGridDB` spatial database, and we specify that we want the following fields to be output for postprocessing: `[displacement, pressure, cauchy_stress, velocity, porosity, isotropic_permeability, water_content]`.

```{code-block} cfg
---
caption: Material parameters and data fields for output for the poroelastic-outerrise-2d examples.
---
[pylithapp.problem]
materials = [slab]
materials.slab = pylith.materials.Poroelasticity

[pylithapp.problem.materials]
slab.bulk_rheology = pylith.materials.IsotropicLinearPoroelasticity

[pylithapp.problem.materials.slab]

description = slab
label_value = 1
use_state_variables = True
db_auxiliary_field = spatialdata.spatialdb.SimpleGridDB 
db_auxiliary_field.description = Test 
db_auxiliary_field.query_type = linear

observers.observer.data_fields = [displacement, pressure, cauchy_stress, velocity, porosity, isotropic_permeability, water_content]
```

For all steps, the left boundary is held fixed, which allows the rest of the boundary to pivot about this anchor point. The top boundary always has a fluid pressure boundary condition, while the displacement boundary condition varies between steps for the top boundary. The bottom and right boundaries remain unconstrained. 

```{code-block} cfg
---
caption: Boundary conditions for the poroelastic-outerrise-2d examples.
---
[pylithapp.problem]
bc = [bndry_west, bndry_top, bndry_top_fluid]

bc.bndry_west = pylith.bc.DirichletTimeDependent
bc.bndry_top = pylith.bc.DirichletTimeDependent
bc.bndry_top_fluid = pylith.bc.DirichletTimeDependent

[pylithapp.problem.bc.bndry_west]
constrained_dof = [0, 1]
label = bndry_west
label_value = 11
field = displacement
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC -x

[pylithapp.problem.bc.bndry_top]
constrained_dof = [0, 1]
label = bndry_top
label_value = 10
field = displacement
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC +y

[pylithapp.problem.bc.bndry_top_fluid]
constrained_dof = [0]
label = bndry_top_fluid
label_value = 14
field = pressure
use_initial = True
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Dirichlet BC +y
db_auxiliary_field.iohandler.filename = simpleDB_files/surface_fluid_pressure.txt
```
