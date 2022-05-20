# Distributor

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.topology.Distributor`
:Journal name: `mesh_distributor`

Distributor of the the mesh among processes.

## Pyre Facilities

* `data_writer`: Data writer for partition information.
  - **current value**: 'datawritervtk', from {default}
  - **configurable as**: datawritervtk, data_writer

## Pyre Properties

* `partitioner`=\<str\>: Name of mesh partitioner.
  - **default value**: 'chaco'
  - **current value**: 'chaco', from {default}
  - **validator**: (in ['chaco', 'metis', 'parmetis', 'simple'])
* `write_partition`=\<bool\>: Write partition information to file.
  - **default value**: False
  - **current value**: False, from {default}

## Example

Example of setting `Distributor` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.mesh_generator.distributor]
partitioner = parmetis
:::

