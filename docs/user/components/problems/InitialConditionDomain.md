# InitialConditionDomain

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.InitialConditionDomain`
:Journal name: `initialconditionsdomain`

Initial conditions for the solution over the entire domain.

Implements `InitialCondition`.

## Pyre Facilities

* `db`: Spatial database with values for initial condition.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db

## Pyre Properties

* `subfields`=\<list\>: Names of solution subfields for initial condition.
  - **default value**: ['displacement']
  - **current value**: ['displacement'], from {default}

## Example

Example of setting `InitialConditionDomain` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Create a single initial condition over the domain.
[pylithapp.problem]
ic = [domain]
ic.domain = pylith.problems.InitialConditionDomain

[pylithapp.problem.ic.domain]
db = spatialdata.spatialdb.SimpleGridDB
db.description = Initial conditions over domain
db.filename = sheardisp_ic.spatialdb
:::

