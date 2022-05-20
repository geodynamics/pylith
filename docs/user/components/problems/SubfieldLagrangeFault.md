# SubfieldLagrangeFault

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.SubfieldLagrangeFault`
:Journal name: `subfieldlagrangefault`

Object for defining attributes of the fault Lagrange multiplier solution subfield.

Implements `SolutionSubfield`.

## Pyre Properties

* `alias`=\<str\>: Name for subfield.
  - **default value**: ''
  - **current value**: 'lagrange_multiplier_fault', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/pyre/inventory/ConfigurableClass.py', line=26, function='__set__'}
  - **validator**: <function validateAlias at 0x124c314c0>
* `basis_order`=\<int\>: Order of basis functions.
  - **default value**: 1
  - **current value**: 1, from {default}
* `cell_basis`=\<str\>: Type of cell basis functions (simplex, tensor, or default). Default is to use type matching cell type.
  - **default value**: 'default'
  - **current value**: 'default', from {default}
  - **validator**: (in ['simplex', 'tensor', 'default'])
* `dimension`=\<int\>: Topological dimension associated with subfield (=-1 will use dimension of domain).
  - **default value**: -1
  - **current value**: -1, from {default}
* `finite_element_space`=\<str\>: Finite-element space (polynomial or point). Point space corresponds to delta functions at quadrature points.
  - **default value**: 'polynomial'
  - **current value**: 'polynomial', from {default}
  - **validator**: (in ['polynomial', 'point'])
* `is_basis_continous`=\<bool\>: Is basis continuous?
  - **default value**: True
  - **current value**: True, from {default}
* `quadrature_order`=\<int\>: Order of numerical quadrature.
  - **default value**: -1
  - **current value**: -1, from {default}

## Example

Example of setting `SubfieldLagrangeFault` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problems.solution.subfields.lagrange_fault]
alias = lagrange_multiplier_fault
basis_order = 1
:::

