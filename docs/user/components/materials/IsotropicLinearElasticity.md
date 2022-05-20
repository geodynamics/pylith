# IsotropicLinearElasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.IsotropicLinearElasticity`
:Journal name: `isotropiclinearelasticity`

Isotropic linear elastic bulk rheology.

Implements `RheologyElasticity`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for physical properties and state variables.
  - **current value**: 'auxiliary_subfields', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/pyre/inventory/ConfigurableClass.py', line=26, function='__set__'}
  - **configurable as**: auxiliary_subfields

## Pyre Properties

* `use_reference_state`=\<bool\>: Use reference stress/strain state.
  - **default value**: False
  - **current value**: False, from {default}

## Example

Example of setting `IsotropicLinearElasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_elastic.rheology]
use_reference_state = False

auxiliary_subfields.shear_modulus.basis_order = 0
auxiliary_subfields.bulk_modulus.basis_order = 0
:::

