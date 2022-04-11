# AbsorbingDampers

% WARNING: Do not edit; this is a generated file!
Full name: `pylith.bc.AbsorbingDampers`

Absorbing dampers boundary condition.

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

## Pyre Properties

* `field`=\<str\>: Solution subfield associated with boundary condition.
  - **default value**: 'displacement'
  - **current value**: 'displacement', from {default}
* `label`=\<str\>: Label identifying boundary.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateLabel at 0x117dfa700>

## Example

Example of setting `AbsorbingDampers` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[bc]
label = boundary_xpos
field = velocity

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.label = Material properties for absorbing boundary
db_auxiliary_field.values = [density, vs, vp]
db_auxiliary_field.data = [2500*kg/m**3, 1.0*km/s, 1.732*km/s]

auxiliary_subfields.density.basis_order = 0
auxiliary_subfields.vp.basis_order = 0
auxiliary_subfields.vs.basis_order = 0            
:::

