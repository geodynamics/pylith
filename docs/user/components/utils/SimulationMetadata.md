# SimulationMetadata

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.utils.SimulationMetadata`
:Journal name: `metadata`

Metadata for simulation.

When using `base` to specify other files with metadata, the other files will append to the `keywords` and `features` lists, whereas other metadata will be overwritten (the same behavior as other Pyre properties).

## Pyre Properties

* `arguments`=\<list\>: Command line arguments for running simulation.
  - **default value**: []
  - **current value**: [], from {default}
* `authors`=\<list\>: Creator(s) of simulation.
  - **default value**: []
  - **current value**: [], from {default}
* `base`=\<list\>: Parameter files with metadata that complement this metadata.
  - **default value**: []
  - **current value**: [], from {default}
* `description`=\<str\>: Description of simulation.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function notEmptyString at 0x124cdf820>
* `features`=\<list\>: PyLith features used in simulation.
  - **default value**: []
  - **current value**: [], from {default}
* `keywords`=\<list\>: Keywords describing simulation.
  - **default value**: []
  - **current value**: [], from {default}
* `pylith_version`=\<list\>: PyLith versions compatible with simulation input files.
  - **default value**: []
  - **current value**: [], from {default}
* `version`=\<str\>: Version number for simulation.
  - **default value**: ''
  - **current value**: '', from {default}

## Example

Example of setting `SimulationMetadata` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.metadata]
base = [pylithapp.cfg]
description = Axial extension using Dirichlet boundary conditions.
keywords = [example, 2D, box, axial extension]
features = [
    Quadrilateral cells,
    pylith.meshio.MeshIOAscii,
    pylith.problems.TimeDependent,
    pylith.materials.Elasticity,
    pylith.materials.IsotropicLinearElasticity,
    spatialdata.spatialdb.UniformDB,
    pylith.meshio.DataWriterHDF5
    ]
authors = [Brad Aagaard]
version = 1.0.0
arguments = [step01_axialdisp.cfg]
pylith_version = [>=3.0, <4.0]
:::

