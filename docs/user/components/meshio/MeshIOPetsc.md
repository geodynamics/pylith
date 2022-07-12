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
* `options_prefix`=\<str\>: Name of PETSc options prefix for this mesh.
  - **default value**: ''
  - **current value**: '', from {default}

