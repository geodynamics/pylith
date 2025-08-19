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

We specify a basis order of 2 for the displacement, the velocity, and the isotropic permeability subfields, and a basis order of 1 for each of the remaining subfields.

```{code-block} cfg
---
caption: Discretization parameters for the 2D outer-rise examples with poroelasticity.
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

```{code-block} cfg
---
caption: Nondimensionalization parameters for the 2D outer-rise examples with poroelasticity.
---
[pylithapp.problem]
scales = spatialdata.units.NondimElasticQuasistatic
scales.length_scale = 100.0*m
scales.relaxation_time = 1*year
scales.shear_modulus = 10.0*GPa
```

```{code-block} cfg
---
caption: Time stepping parameters for the 2D outer-rise examples with poroelasticity.
---
[pylithapp.timedependent]
start_time = -6e3*year
initial_dt = 6e3*year
end_time = 300e3*year
```

We set the material parameters that are common across the three steps in the `pylithapp.cfg` file to avoid repeating parameters.
We only have a single material with spatially varying material properties so we use a `SimpleGridDB` spatial database.
We use differ `SimpleGridDB` files to initialize the material properties for the three steps, so we specify the filename in the parameter file for each step.

```{code-block} cfg
---
caption: Material parameters common to all steps in the 2D outer-rise examples with poroelasticity.
---
[pylithapp.problem]
materials = [slab]
materials.slab = pylith.materials.Poroelasticity

[pylithapp.problem.materials]
slab.bulk_rheology = pylith.materials.IsotropicLinearPoroelasticity

[pylithapp.problem.materials.slab]
label_value = 1
use_state_variables = True
db_auxiliary_field = spatialdata.spatialdb.SimpleGridDB 
db_auxiliary_field.description = Spatial database for material properties and state variables
db_auxiliary_field.query_type = linear

observers.observer.data_fields = [displacement, pressure, cauchy_stress, velocity, porosity, isotropic_permeability, water_content]
```

For all steps, the left boundary is held fixed, which allows the rest of the boundary to pivot about this anchor point.
The top boundary always has a fluid pressure boundary condition, while the displacement boundary condition varies between steps for the top boundary.
The bottom and right boundaries remain unconstrained.

```{code-block} cfg
---
caption: Boundary conditions for the 2D outer-rise examples with poroelasticity.
---
[pylithapp.problem]
bc = [bc_xneg, bc_ypos, bc_ypos_fluid]

bc.bc_west = pylith.bc.DirichletTimeDependent
bc.bc_ypos = pylith.bc.DirichletTimeDependent
bc.bc_ypos_fluid = pylith.bc.DirichletTimeDependent

[pylithapp.problem.bc.bc_xneg]
constrained_dof = [0, 1]
label = boundary_xneg
label_value = 11
field = displacement
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC on -x for displacement

[pylithapp.problem.bc.bc_ypos]
constrained_dof = [0, 1]
label = boundary_ypos
label_value = 10
field = displacement
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC on +y for displacement

[pylithapp.problem.bc.bc_ypos_fluid]
constrained_dof = [0]
label = boundary_ypos
label_value = 10
field = pressure
use_initial = True
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Dirichlet BC on +y for fluid pressure
db_auxiliary_field.iohandler.filename = simpleDB_files/surface_fluid_pressure.txt
```

## Generate spatial databases

:::{important}
The spatial database files are not provided because some of them are large.
You generate them using the provided Python scripts _before_ running the simulations.
:::

```{code-block} console
---
caption: Generate the spatial database files.
---
# Generate spatial databases for the Dirichlet BC on the top surface.
./generate_spatialdb_ypos.py

# Generate spatial databases for the material properties.
./generate_spatialdb_matprops.py
```
