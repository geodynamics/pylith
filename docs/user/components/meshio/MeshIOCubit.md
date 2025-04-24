# MeshIOCubit

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.MeshIOCubit`
:Journal name: `meshiocubit`

Reader for finite-element meshes from Exodus II files (usually from Cubit).

:::{warning}
The coordinate system associated with the mesh must be a Cartesian coordinate system, such as a generic Cartesian coordinate system or a geographic projection.
:::

Implements `MeshIOObj`.

## Pyre Facilities

* `coordsys`: Coordinate system associated with mesh.
  - **current value**: 'cscart', from {default}
  - **configurable as**: cscart, coordsys

## Pyre Properties

* `filename`=\<str\>: Name of Cubit Exodus file.
  - **default value**: 'mesh.exo'
  - **current value**: 'mesh.exo', from {default}

## Example

Example of setting `MeshIOCubit` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.mesh_generator.reader]
filename = mesh_quad.exo
coordsys.space_dim = 2
:::

