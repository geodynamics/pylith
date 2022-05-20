# InitialConditionPatch

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.InitialConditionPatch`
:Journal name: `initialconditionspatch`

Initial conditions over a portion of the domain (patch).

Implements `InitialCondition`.

## Pyre Facilities

* `db`: Spatial database with values for initial condition.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db

## Pyre Properties

* `label`=\<str\>: Name of label for patch.
  - **default value**: 'material-id'
  - **current value**: 'material-id', from {default}
* `label_value`=\<int\>: Value of label associated with initial condition patch, usually the material label value.
  - **default value**: 1
  - **current value**: 1, from {default}
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
label_value = 1
db = spatialdata.spatialdb.SimpleGridDB
db.description = Initial conditions over material 1
db.filename = shearmat1_ic.spatialdb

[pylithapp.problem.ic.mat2]
label_value = 2
db = spatialdata.spatialdb.SimpleGridDB
db.description = Initial conditions over material 2
db.filename = shearmat2_ic.spatialdb
:::

