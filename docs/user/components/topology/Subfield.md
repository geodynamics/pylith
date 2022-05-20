# Subfield

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.topology.Subfield`
:Journal name: `subfield`

Python object for defining discretization of a subfield.

## Pyre Properties

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

