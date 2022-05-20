# DirichletTimeDependent

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.bc.DirichletTimeDependent`
:Journal name: `dirichlettimedependent`

Dirichlet (prescribed values) time-dependent boundary condition.

This boundary condition sets values of a single solution subfield on a boundary.
To set multiple solution subfields on a boundary, use multiple Dirichlet boundary conditions.

:::{seealso}
See [`AuxSubfieldsTimeDependent` Component](AuxSubfieldsTimeDependent.md) for the functional form of the time depenence.
:::

Implements `BoundaryCondition`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'auxiliary_subfields', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/pyre/inventory/ConfigurableClass.py', line=26, function='__set__'}
  - **configurable as**: auxiliary_subfields
* `db_auxiliary_field`: Database for physical property parameters.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db_auxiliary_field
* `derived_subfields`: Discretization of derived subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, derived_subfields
* `observers`: Observers (e.g., output).
  - **current value**: 'singlephysicsobserver', from {default}
  - **configurable as**: singlephysicsobserver, observers
* `time_history`: Time history with normalized amplitude.
  - **current value**: 'nullcomponent', from {default}
  - **configurable as**: nullcomponent, time_history

## Pyre Properties

* `constrained_dof`=\<array\>: Array of constrained degrees of freedom (0=1st DOF, 1=2nd DOF, etc).
  - **default value**: []
  - **current value**: [], from {default}
* `field`=\<str\>: Solution subfield associated with boundary condition.
  - **default value**: 'displacement'
  - **current value**: 'displacement', from {default}
* `label`=\<str\>: Name of label identifying boundary.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateLabel at 0x124bbc4c0>
* `label_value`=\<int\>: Value of label identifying boundary (tag of physical group in Gmsh files).
  - **default value**: 1
  - **current value**: 1, from {default}
* `use_initial`=\<bool\>: Use initial term in time-dependent expression.
  - **default value**: True
  - **current value**: True, from {default}
* `use_rate`=\<bool\>: Use rate term in time-dependent expression.
  - **default value**: False
  - **current value**: False, from {default}
* `use_time_history`=\<bool\>: Use time history term in time-dependent expression.
  - **default value**: False
  - **current value**: False, from {default}

## Example

Example of setting `DirichletTimeDependent` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Dirichlet (prescribed displacements) boundary condition constraining the x and y degrees of freedom on the +y boundary.
[pylithapp.problem.bc.bc_ypos]
constrained_dof = [0, 1]
label = boundary_ypos
field = displacement

use_initial = False
use_time_history = True
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Displacement Dirichlet BC +y boundary
db_auxiliary_field.values = [time_history_amplitude_x, time_history_amplitude_y, time_history_start_time]
db_auxiliary_field.data = [1.0*m, 0.0*m, 0.0]

time_history = spatialdata.spatialdb.TimeHistory
time_history.description = Impulse time history
time_history.filename = impulse.timedb
:::

