# MeshIOPetsc

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.MeshIOPetsc`
:Journal name: `meshiopetsc`

Python object for a variety of reading/writing finite-element meshes using PETSc.
Currently, the primary use of this object is to import meshes from Gmsh.

:::{warning}
The coordinate system associated with the mesh must be a Cartesian coordinate system, such as a generic Cartesian coordinate system or a geographic projection.
:::

Implements `MeshIOObj`.

## Pyre Facilities

* `coordsys`: Coordinate system associated with mesh.
  - **current value**: 'cscart', from {default}
  - **configurable as**: cscart, coordsys

## Pyre Properties

* `filename`=\<str\>: Name of mesh file for reading with PETSc.
  - **default value**: ''
  - **current value**: '', from {default}
* `gmsh_mark_vertices`=\<bool\>: Gmsh file marks faces, edges, and vertices rather than just faces.
  - **default value**: False
  - **current value**: False, from {default}
* `options_prefix`=\<str\>: Name of PETSc options prefix for this mesh.
  - **default value**: ''
  - **current value**: '', from {default}

## Example

Example of setting `MeshIOPetsc` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.mesh_generator.reader]
filename = mesh_quad.msh
gmsh_mark_vertices = False
coordsys.space_dim = 2
:::

