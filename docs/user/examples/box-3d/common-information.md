# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`.
This limits duplicate information in the parameter files for each simulation.

## Metadata, Mesh, and Output

The `pylithapp.metadata` section specifies metadata common to all simulations in the directory.
We control the verbosity of the output written to stdout using `journal.info`.
We use `MeshIOPetsc` to import finite-element mesh from Gmsh.

```{code-block} cfg
---
caption: Parameters for importing the finite-element mesh from Gmsh.
---
[pylithapp.mesh_generator]
reader = pylith.meshio.MeshIOPetsc

# Set the filename and the dimension of the default Cartesian coordinate system.
reader.filename = mesh_hex.msh
reader.coordsys.space_dim = 3
```

In addition to output of the solution over the domain (default), we also output the solution on the ground surface (+z boundary).

```{code-block} cfg
---
caption: Parameters for output over the domain and ground surface (+z boundary).
---
[pylithapp.problem]
solution_observers = [domain, ground_surface]
solution_observers.ground_surface = pylith.meshio.OutputSolnBoundary

[pylithapp.problem.solution_observers.ground_surface]
# The `label` and `label_value` match the name and tag in the Gmsh Python script.
label = boundary_zpos
label_value = 15
```

## Physics

All of the simulations in this example suite solve the elasticity equation without inertia,
%
\begin{gather}
\vec{s} = \left(\begin{array}{c} \vec{u} \end{array}\right)^T \\
\boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u}) = \vec{0}
\end{gather}
%
The default PyLith settings target solving the quasi-static elasticity equation, so we can use the default `TimeDependent` problem and solution field with a single displacement subfield of basis order 1.

In this suite of examples we have a single material, and we must specify a description for the material and the label value.
In this case the label value corresponds to the tag of the physical group in the Gmsh file.
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
