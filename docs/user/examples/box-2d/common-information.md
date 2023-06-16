# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`.
This limits duplicate information in the parameter files for each simulation.

## Metadata, Mesh, and Output

The `pylithapp.metadata` section specifies metadata common to all simulations in the directory.
We control the verbosity of the output written to stdout using `journal.info`.
We set the parameters for importing the finite-element mesh in `pylithapp.mesh_generator`. 

## Physics

All of the simulations in this example suite solve the elasticity equation without inertia,
%
\begin{gather}
\vec{s} = \left(\begin{array}{c} \vec{u} \end{array}\right)^T \\
\boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u}) = \vec{0}
\end{gather}
%
The default PyLith settings target solving the quasi-static elasticity equation, so we can use the default `TimeDependent` problem and solution field with a single displacement subfield of basis order 1.

```{code-block} cfg
---
caption: Default PyLith problem and solution field settings appropriate for quasi-static elasticity. We do not include these in `pylithapp.cfg`.
---
[pylithapp]
# Default problem and solution field.
problem = pylith.problems.TimeDependent
problem.solution = pylith.problems.SolnDisp
        
[pylithapp.problem.solution.subfields]
# The default basis order for solution subfields is 1 (linear variation within cells).
displacement.basis_order = 1

[pylithapp.problem]
# Default output for the solution is over the domain.
solution_observers = [domain]
```

The default material and bulk rheology also target solving the quasi-static elasticity equation.
The default bulk rheology is isotropic linear elasticity.

```{code-block} cfg
---
caption: Default PyLith material and bulk rheology settings appropriate for quasi-static elasticity. We do not include these in `pylithapp.cfg`.
---
[pylithapp.problem]
materials = [elastic]
materials.elastic = pylith.materials.Elasticity

[pylithapp.problem.materials.elastic]
elastic.bulk_rheology = pylith.materials.IsotropicLinearElasticity
```

In this suite of examples we have a single material, and we must specify a description for the material and the label value.
In this case the label value corresponds to the `material-id` in the finite-element mesh file.
The physical properties for each material are specified in a spatial database.
The material properties are uniform, so we use a single `UniformDB` to specify the material properties in `pylithapp.cfg`.
When we have a spatial variation in the material properties, we specify them in a spatial database file.
With uniform properties we use a basis order of 0 for the auxiliary subfields.

```{code-block} cfg
---
caption: Material settings in `pylithapp.cfg`.
---
[pylithapp.problem]
materials = [elastic]

[pylithapp.problem.materials.elastic]
description = Elastic material
label_value = 0

# We will use uniform material properties, so we use the UniformDB
# spatial database.
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Elastic properties
db_auxiliary_field.values = [density, vs, vp]
db_auxiliary_field.data = [2500*kg/m**3, 3.0*km/s, 5.2915026*km/s]

# Set the discretization of the material auxiliary fields (properties).
# We have uniform material properties, so we can use a basis order of 0.
auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0
```
