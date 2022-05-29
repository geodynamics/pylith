# PetscDefaults

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.utils.PetscDefaults`
:Journal name: `petscdefaults`

Flags controlling use of default PETSc settings.
No user-specified settings will be overwritten.

## Pyre Properties

* `monitors`=\<bool\>: Use default solver monitors.
  - **default value**: True
  - **current value**: True, from {default}
* `parallel`=\<bool\>: Use solver settings normally used when running in parallel.
  - **default value**: False
  - **current value**: False, from {default}
* `solver`=\<bool\>: Use default solver settings based on governing equations.
  - **default value**: True
  - **current value**: True, from {default}
* `testing`=\<bool\>: Use default PETSc testingging options.
  - **default value**: False
  - **current value**: False, from {default}

## Example

Example of setting `PetscDefaults` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.petsc_defaults]
solver = True
monitors = True
parallel = False
testing = False
:::

