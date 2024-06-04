# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`, which contains numerous comments, so we only summarize the parameters here.

## Metadata, Mesh, and Output

The `pylithapp.metadata` section specifies metadata common to all simulations in the directory.
We control the verbosity of the output written to stdout using `journal.info`.
We set the parameters for importing the finite-element mesh in `pylithapp.mesh_generator`.

:::{important}
We specify the geographic projection that we used when we constructed the finite-element mesh.
This georeferences the model and enables us to specify points in other geographic coordinates system in spatial databases, which may be more convenient.
:::

```{code-block} cfg
---
caption: Parameters for the mesh input file and the model coordinate system.
---
[pylithapp.mesh_generator]
reader = pylith.meshio.MeshIOPetsc
reader.filename = mesh_tri.msh

# Model coordinate system is UTM zone 11
reader.coordsys = spatialdata.geocoords.CSGeo
reader.coordsys.crs_string = EPSG:32611
reader.coordsys.space_dim = 3
```

## Physics

These quasi-static simulations solve the elasticity equation and include a fault, so we have a solution field with both displacement and Lagrange multiplier subfields.
%
\begin{gather}
\vec{s} = \left(\begin{array}{c} \vec{u} \quad \vec{\lambda} \end{array}\right)^T \\
\boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u}) = \vec{0}
\end{gather}
%
We use the default `TimeDependent` problem and solution field with a single displacement subfield of basis order 1.

```{code-block} cfg
---
caption: Parameters for the solution and output of the solution on the +z boundary (ground surface).
---
[pylithapp.problem]
solution = pylith.problems.SolnDispLagrange

[pylithapp.problem]
# Output the solution over the domain and on the +z (groundsurf).
solution_observers = [domain, groundsurf]
solution_observers.groundsurf = pylith.meshio.OutputSolnBoundary

# `label` and `label_value` correspond to the name and tag of the 
# physical groups marking the boundaries in the Gmsh script.
[pylithapp.problem.solution_observers.groundsurf]
label = boundary_top
label_value = 15
```

We use the same material properties in all of the simulations in this directory, so we specify them in `pylithapp.cfg` to avoid repeating the information in the file with parameters for each simulation.
We have one material with uniform material properties.

```{code-block} cfg
---
caption: Material parameters common to all simulations in this directory.
---
[pylithapp.problem]
materials = [elastic]

[pylithapp.problem.materials.elastic]
# We use the default material (elasticity) and rheology
# (isotropic, linearly elastic).
#
# `label_value` must match the tag for the physical group in the Gmsh Python script.
description = Elastic material
label_value = 1

# The properties are uniform within the material, so we use a UniformDB.
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Elastic properties xneg
db_auxiliary_field.values = [density, vs, vp]
db_auxiliary_field.data = [2500.0*kg/m**3, 3.00*km/s, 5.29*km/s]

# The properties are uniform, so we use a basis order of 0, corresponding
# to uniform properties within a cell.
auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0

# We discretize the displacement field with a basis order of 1
# so the stress and strain computed from the displacement field
# will have an accuracy of one order lower. Consequently, we use
# a basis order of 0.
derived_subfields.cauchy_strain.basis_order = 0
derived_subfields.cauchy_stress.basis_order = 0
```

Similarly, we set the fault parameters common to all simulations in the directory.
For these simulations, we have three faults that intersect at a common point.

```{code-block} cfg
---
caption: General fault parameters common to all simulations in this directory.
---
[pylithapp.problem]
interfaces = [main_fault, west_branch, east_branch]

## main fault
[pylithapp.problem.interfaces.main_fault]
# The `label` and `label_value` correspond to the name and tag of the
# physical group in the Gmsh Python script.
label = fault_main
label_value = 20

edge = fault_main_edges
edge_value = 30

# Output `slip` and `traction_changes` on the fault.
observers.observer.data_fields = [slip, traction_change]


[pylithapp.problem.interfaces.west_branch]
# The `label` and `label_value` correspond to the name and tag of the
# physical group in the Gmsh Python script.
label =  fault_west
label_value = 21

edge = fault_west_edges 
edge_value = 31

# Output `slip` and `traction_changes` on the fault.
observers.observer.data_fields = [slip, traction_change]


[pylithapp.problem.interfaces.east_branch]
# The `label` and `label_value` correspond to the name and tag of the
## physical group in the Gmsh Python script.
label = fault_east
label_value = 22

edge = fault_east_edges
edge_value = 32

# Output `slip` and `traction_changes` on the fault.
observers.observer.data_fields = [slip, traction_change]
```

All simulations in this directory use Dirichlet boundary conditions on the lateral and bottom boundaries of the domain with zero displacements in which we constrain the displacement component perpendicular to the boundary.

```{code-block} cfg
---
caption: Dirichlet boundary conditions common to all simulations in this directory.
---
[pylithapp.problem]
bc = [bc_south, bc_east, bc_north, bc_west, bc_bottom]
bc.bc_south = pylith.bc.DirichletTimeDependent
bc.bc_east = pylith.bc.DirichletTimeDependent
bc.bc_north = pylith.bc.DirichletTimeDependent
bc.bc_south = pylith.bc.DirichletTimeDependent
bc.bc_bottom = pylith.bc.DirichletTimeDependent

# The `label` and `label_value` correspond to the name and tag of the
# physical group in the Gmsh Python script.
#
# We constrain the displacement component normal to each boundary.
# We use the specialized `ZeroDB` to specify zero values for the Dirichlet
# BC. We will override this parameter in some of the .cfg files to specify
# nonzero values.
[pylithapp.problem.bc.bc_south]
label = boundary_south
label_value = 10
constrained_dof = [1]
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC on south boundary

[pylithapp.problem.bc.bc_east]
label = boundary_east
label_value = 11
constrained_dof = [0]
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC on east boundary

[pylithapp.problem.bc.bc_north]
label = boundary_north
label_value = 12
constrained_dof = [1]
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC on north boundary

[pylithapp.problem.bc.bc_west]
label = boundary_west
label_value = 13
constrained_dof = [0]
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC on west boundary

[pylithapp.problem.bc.bc_bottom]
label = boundary_bottom
label_value = 14
constrained_dof = [2]
db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet BC on bottom boundary
```
