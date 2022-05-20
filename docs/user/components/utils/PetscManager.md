# PetscManager

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.utils.PetscManager`
:Journal name: `petsc`

Manage PETSc options.

## Example

Example of setting `PetscManager` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.petsc]
ts_monitor = true
ksp_monitor = true
ksp_converged_reason = true
snes_monitor = true
snes_converged_reason = true
snes_linesearch_monitor = true
:::

