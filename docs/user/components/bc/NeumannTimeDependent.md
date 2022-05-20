# NeumannTimeDependent

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.bc.NeumannTimeDependent`
:Journal name: `neumanntimedependent`

Neumann time-dependent boundary condition. Implements `BoundaryCondition`.

This boundary condition applies a Neumann boundary condition for a single solution subfield on a boundary.
To apply Neumann boundary conditions for multiple solution subfields on a boundary, use multiple Neumann boundary conditions.

:::{important}
The components are specified in the local normal-tangential coordinate system for the boundary. Ambiguities in specifying the shear (tangential) tractions in 3D problems are resolved using the `ref_dir_1` and `ref_dir_2` properties.
The first tangential direction is $ec{z} 	imes ec{r}_1$ unless these are colinear, then $ec{r}_2$ (`ref_dir_2`) is used.
The second tangential direction is $ec{n} 	imes ec{t}_1$.
:::

:::{seealso}
See [`AuxSubfieldsTimeDependent` Component](AuxSubfieldsTimeDependent.md) for the functional form of the time depenence.
:::

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
* `time_history`: Time history with normalized amplitude as a function of time.
  - **current value**: 'nullcomponent', from {default}
  - **configurable as**: nullcomponent, time_history

## Pyre Properties

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
* `ref_dir_1`=\<list\>: First choice for reference direction to discriminate among tangential directions in 3D.
  - **default value**: [0.0, 0.0, 1.0]
  - **current value**: [0.0, 0.0, 1.0], from {default}
  - **validator**: <function validateDir at 0x124bbc9d0>
* `ref_dir_2`=\<list\>: Second choice for reference direction to discriminate among tangential directions in 3D.
  - **default value**: [0.0, 1.0, 0.0]
  - **current value**: [0.0, 1.0, 0.0], from {default}
  - **validator**: <function validateDir at 0x124bbc9d0>
* `scale_name`=\<str\>: Type of scale for nondimensionalizing Neumann boundary condition ('pressure' for elasticity).
  - **default value**: 'pressure'
  - **current value**: 'pressure', from {default}
  - **validator**: (in ['length', 'time', 'pressure', 'density', 'velocity'])
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

Example of setting `NeumannTimeDependent` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Neumann (traction) boundary condition in 2D on -y boundary.
[pylithapp.problem.bc.bc_yneg]
label = boundary_yneg
field = displacement
scale_name = pressure

use_initial = False
use_time_history = True
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Displacement Neumann BC +y boundary
db_auxiliary_field.values = [time_history_amplitude_tangential, time_history_amplitude_normal, time_history_start_time]
db_auxiliary_field.data = [2.0*MPa, -1.0*MPa, 0.0]

time_history = spatialdata.spatialdb.TimeHistory
time_history.description = Impulse time history
time_history.filename = impulse.timedb
:::

