# Problem

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.problems.Problem`
:Journal name: `problem`

Abstract base class for a problem.

## Pyre Facilities

* `bc`: Boundary conditions.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, bc
* `defaults`: Default options for problem.
  - **current value**: 'problem_defaults', from {default}
  - **configurable as**: problem_defaults, defaults
* `gravity_field`: Database used for gravity field.
  - **current value**: 'nullcomponent', from {default}
  - **configurable as**: nullcomponent, gravity_field
* `interfaces`: Interior surfaces with constraints or constitutive models.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, interfaces
* `materials`: Materials in problem.
  - **current value**: 'homogeneous', from {default}
  - **configurable as**: homogeneous, materials
* `normalizer`: Nondimensionalizer for problem.
  - **current value**: 'nondimelasticquasistatic', from {default}
  - **configurable as**: nondimelasticquasistatic, normalizer
* `petsc_defaults`: Flags controlling which default PETSc options to use.
  - **current value**: 'petscdefaults', from {default}
  - **configurable as**: petscdefaults, petsc_defaults
* `solution`: Solution field for problem.
  - **current value**: 'solution', from {default}
  - **configurable as**: solution
* `solution_observers`: Observers (e.g., output) for solution.
  - **current value**: 'singlesolnobserver', from {default}
  - **configurable as**: singlesolnobserver, solution_observers

## Pyre Properties

* `formulation`=\<str\>: Formulation for equations.
  - **default value**: 'quasistatic'
  - **current value**: 'quasistatic', from {default}
  - **validator**: (in ['quasistatic', 'dynamic', 'dynamic_imex'])
* `solver`=\<str\>: Type of solver to use ['linear', 'nonlinear'].
  - **default value**: 'linear'
  - **current value**: 'linear', from {default}
  - **validator**: (in ['linear', 'nonlinear'])

