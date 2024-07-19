# Distributor

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.topology.Distributor`
:Journal name: `mesh_distributor`

Distributor of the the mesh among processes.

## Pyre Facilities

* `data_writer`: Data writer for partition information.
  - **current value**: 'datawriterhdf5', from {default}
  - **configurable as**: datawriterhdf5, data_writer

## Pyre Properties

* `partitioner`=\<str\>: Name of mesh partitioner (PETSc must be built with partitioner).
  - **default value**: 'parmetis'
  - **current value**: 'parmetis', from {default}
  - **validator**: (in ['parmetis', 'chaco', 'simple'])
* `use_edge_weighting`=\<bool\>: Use edge weighting (parmetis only).
  - **default value**: True
  - **current value**: True, from {default}
* `write_partition`=\<bool\>: Write partition information to file.
  - **default value**: False
  - **current value**: False, from {default}

## Example

Example of setting `Distributor` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.mesh_generator.distributor]
partitioner = parmetis
:::

