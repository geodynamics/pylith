# InitialConditionPatch

% WARNING: Do not edit; this is a generated file!
Full name: `pylith.problems.InitialConditionPatch`

Initial conditions over a portion of the domain (patch).

Implements `InitialCondition`.

## Pyre Facilities

* `db`: Spatial database with values for initial condition.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db

## Pyre Properties

* `id`=\<int\>: Material id associated with patch.
  - **default value**: 0
  - **current value**: 0, from {default}
* `subfields`=\<list\>: Names of solution subfields for initial condition.
  - **default value**: ['displacement']
  - **current value**: ['displacement'], from {default}

## Example

Example of setting `InitialConditionPatch` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Create separate initial conditions for two materials.
# This is often useful if the materials have different properties.
[pylithapp.problem]
ic = [mat1, mat2]
ic.mat1 = pylith.problems.InitialConditionPatch
ic.mat2 = pylith.problems.InitialConditionPatch

[pylithapp.problem.ic.mat1]
id = 1
db = spatialdata.spatialdb.SimpleGridDB
db.label = Initial conditions over material 1
db.filename = shearmat1_ic.spatialdb

[pylithapp.problem.ic.mat2]
id = 2
db = spatialdata.spatialdb.SimpleGridDB
db.label = Initial conditions over material 2
db.filename = shearmat2_ic.spatialdb
:::

