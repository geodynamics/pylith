# ZeroDB

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.bc.ZeroDB`
:Journal name: `zerodb`

Special case of a `UniformDB` spatial database with uniform zero initial amplitude values for degrees of freedom.

Implements `SpatialDB`.

## Pyre Properties

* `data`=\<list\>: Values in spatial database.
  - **default value**: []
  - **current value**: [], from {default}
* `description`=\<str\>: Description for database.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateDescription at 0x124922820>
* `label`=\<str\>: Label for ZeroDB spatial database.
  - **default value**: 'Zero initial amplitude spatial database.'
  - **current value**: 'Zero initial amplitude spatial database.', from {default}
* `values`=\<list\>: Names of values in spatial database.
  - **default value**: []
  - **current value**: [], from {default}

## Example

Example of setting `ZeroDB` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
# Zero displacement boundary condition constraining the y degree of freedom on the -y boundary.
[pylithapp.problem.bc.bc_yneg]
constrained_dof = [1]
label = boundary_yneg
field = displacement

db_auxiliary_field = pylith.bc.ZeroDB
db_auxiliary_field.description = Dirichlet displacement boundary condition on the -y boundary
:::

