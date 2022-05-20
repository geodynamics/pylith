# IsotropicLinearPoroelasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.materials.IsotropicLinearPoroelasticity`
:Journal name: `isotropiclinearporoelasticity`

Isotropic linear incompressible elastic bulk rheology.

Implements `RheologyIncompressibleElasticity`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for physical properties and state variables.
  - **current value**: 'auxiliary_subfields', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/pyre/inventory/ConfigurableClass.py', line=26, function='__set__'}
  - **configurable as**: auxiliary_subfields

## Pyre Properties

* `use_reference_state`=\<bool\>: Use reference stress/strain state.
  - **default value**: False
  - **current value**: False, from {default}
* `use_tensor_permeability`=\<bool\>: Use tensor permeability.
  - **default value**: False
  - **current value**: False, from {default}

## Example

Example of setting `IsotropicLinearPoroelasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.materials.mat_incompelastic.rheology]
use_reference_state = False

auxiliary_subfields.shear_modulus.basis_order = 0
auxiliary_subfields.bulk_modulus.basis_order = 0
:::

