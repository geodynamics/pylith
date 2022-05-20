# ProblemDefaults

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.ProblemDefaults`
:Journal name: `problem_defaults`

Default options for a problem.
Specifying defaults at the problem level (here) will override defaults for individual components.
Non-default values specified for individual components will override the problem defaults (specified here).

## Pyre Properties

* `name`=\<str\>: Name for the problem (used with output_directory for default output filenames).
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateName at 0x124c42280>
* `output_basis_order`=\<int\>: Default basis order for output.
  - **default value**: 1
  - **current value**: 1, from {default}
  - **validator**: (in [0, 1])
* `output_directory`=\<str\>: Directory for output.
  - **default value**: 'output'
  - **current value**: 'output', from {default}
* `quadrature_order`=\<int\>: Finite-element quadrature order.
  - **default value**: 1
  - **current value**: 1, from {default}
  - **validator**: (greater than 0)

## Example

Example of setting `ProblemDefaults` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.defaults]
output_directory = output
name = step01
quadrature_order = 1
output_basis_order = 0
:::

