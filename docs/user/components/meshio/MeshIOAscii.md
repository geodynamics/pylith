# MeshIOAscii

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.MeshIOAscii`
:Journal name: `meshioascii`

Reader for finite-element meshes using a simple ASCII format.

:::{warning}
The coordinate system associated with the mesh must be a Cartesian coordinate system, such as a generic Cartesian coordinate system or a geographic projection.
:::

Implements `MeshIOObj`.

## Pyre Facilities

* `coordsys`: Coordinate system associated with mesh.
  - **current value**: 'cscart', from {default}
  - **configurable as**: cscart, coordsys

## Pyre Properties

* `filename`=\<str\>: Name of mesh file
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateFilename at 0x11f289670>

## Example

Example of setting `MeshIOAscii` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.mesh_generator.reader]
filename = mesh_quad.txt
coordsys.space_dim = 2
:::

