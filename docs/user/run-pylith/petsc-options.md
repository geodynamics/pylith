(sec-user-run-pylith-petsc-options)=
# PETSc Options

PyLith relies on PETSc for finite-element data structures, linear and nonlinear solvers, and time-stepping algorithms.
PETSc has its own object-oriented interface for specifying runtime options.
PyLith provides two mechanisms for passing options to PETSc:

1. Default values that are controlled by a few flags in the PyLith `Problem` `petsc_defaults` facility, and
2. PyLith `petsc` facility which passes options directly to PETSc.

When using the default values, PyLith selects solver and preconditioner options based on the governing equations.
In most cases these default values will give good performance and you do not need to specify any PETSc options.
The user-specified values always take precedence over the default values.

## Default PETSc Options

:::{note}
New in v3.0.0
:::

We separate the defaults into a few categories to make it easy to select desired options.

:solver: Options for the preconditioner and solver;
:parallel: Options used when running in parallel (can be used in serial as well);
:monitors: Options for basic monitoring of the solver; and
:testing: Options used in testing.

:::{tip}
You can see which options PyLith sets using the `petscoptions` Pyre Journal.
Either use the `--journal.info.petscoptions` command line argument or in your `.cfg` file include

```{code-blocK} cfg
[journal.info]
petscoptions = 1
```

:::

:::{seealso}
See [`PetscDefaults` Component](../components/utils/PetscDefaults.md) for more information about the the Pyre interface for specifying default PETSc options.
:::

### Solver Options

The solver options are enabled by default.
PyLith selects options based on the governing equation, formulation, presence of a fault, and whether the simulation is running in parallel.
In most cases the options used when running in parallel give comparable or better performance than those used when running serial; consequently, you may want to use them when solving moderate to larger problems in serial.
Additionally, PyLith specifies general options related to the solver tolerances and triggering errors if the linear or nonlinear solver fails to converge.
The different sets of defaults are detailed in the following code blocks.

:::{warning}
The split fields and algebraic multigrid preconditioning currently fails in problems with a nonzero null space.
This most often occurs when a problem contains multiple faults that extend through the entire domain and create subdomains without any Dirichlet boundary conditions.
The workaround is to use the `ilu` preconditioner.
However, it only works in serial.
An alternative is to use the `asm` preconditioner (Additive Schwarz) which works in parallel and serial.
:::

```{code-block} cfg
---
caption: PETSc options set whenever the solver defaults are enabled.
---
[pylithapp.petsc]
ksp_rtol = 1.0e-12
ksp_atol = 1.0e-12
ksp_error_if_not_converged = true

snes_rtol = 1.0e-12
snes_atol = 1.0e-9
snes_error_if_not_converged = true
```

#### Quasistatic Elasticity

```{code-block} cfg
---
caption: PETSc options used for quasistatic elasticity in serial without a fault.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = lu
```

```{code-block} cfg
---
caption: PETSc options used for quasistatic elasticity in parallel without a fault.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = gamg
mg_levels_pc_type = sor
mg_levels_ksp_type = richardson
```

The Lagrange multiplier corresponding to the tractions on the fault introduces a saddle point in the system of equations, so we use a Schur complement approach.

```{code-block} cfg
---
caption: PETSc options used for quasistatic elasticity in serial with a fault.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_use_amat = true
pc_fieldsplit_type = schur

pc_fieldsplit_schur_factorization_type = lower
pc_fieldsplit_schur_precondition = selfp
pc_fieldsplit_schur_scale = 1.0

fieldsplit_displacement_ksp_type = preonly
fieldsplit_displacement_pc_type = ilu

fieldsplit_lagrange_multiplier_fault_ksp_type = preonly
fieldsplit_lagrange_multiplier_fault_pc_type = lu
```

```{code-block} cfg
---
caption: PETSc options used for quasistatic elasticity in parallel with a fault.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_use_amat = true
pc_fieldsplit_type = schur

pc_fieldsplit_schur_factorization_type = lower
pc_fieldsplit_schur_precondition = selfp
pc_fieldsplit_schur_scale = 1.0

fieldsplit_displacement_ksp_type = preonly
fieldsplit_displacement_pc_type = gamg
fieldsplit_displacement_mg_levels_pc_type = sor
fieldsplit_displacement_mg_levels_ksp_type = richardson

fieldsplit_lagrange_multiplier_fault_ksp_type = preonly
fieldsplit_lagrange_multiplier_fault_pc_type = gamg
fieldsplit_lagrange_multiplier_fault_mg_levels_pc_type = sor
fieldsplit_lagrange_multiplier_fault_mg_levels_ksp_type = richardson
```

```{code-block} cfg
---
caption: Alternative options for quasistatic elasticity with a fault that often provide similar performance.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_use_amat = true
pc_fieldsplit_type = schur

pc_fieldsplit_schur_factorization_type = full
pc_fieldsplit_schur_precondition = selfp
pc_fieldsplit_schur_scale = 1.0

fieldsplit_displacement_ksp_type = preonly
fieldsplit_displacement_pc_type = gamg
fieldsplit_displacement_mg_levels_pc_type = sor
fieldsplit_displacement_mg_levels_ksp_type = richardson

fieldsplit_lagrange_multiplier_fault_ksp_type = preonly
fieldsplit_lagrange_multiplier_fault_pc_type = gamg
fieldsplit_lagrange_multiplier_fault_mg_levels_pc_type = sor
fieldsplit_lagrange_multiplier_fault_mg_levels_ksp_type = richardson
```

#### Quasistatic Incompressible Elasticity

The pressure field introduces a saddle point in the system of equations, so we use a Schur complement approach.

```{code-block} cfg
---
caption: PETSc options used for quasistatic incompressible elasticity in serial.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_fieldsplit_type = schur

pc_fieldsplit_schur_factorization_type = full
pc_fieldsplit_schur_precondition = full

fieldsplit_displacement_pc_type = lu
fieldsplit_pressure_pc_type = lu
```

```{code-block} cfg
---
caption: PETSc options used for quasistatic incompressible elasticity in parallel.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_fieldsplit_type = schur

pc_fieldsplit_schur_factorization_type = full
pc_fieldsplit_schur_precondition = full

fieldsplit_displacement_pc_type = gamg
fieldsplit_displacement_mg_levels_pc_type = sor
fieldsplit_displacement_mg_levels_ksp_type = richardson
fieldsplit_pressure_pc_type = gamg
fieldsplit_pressure_mg_levels_pc_type = sor
fieldsplit_pressure_mg_levels_ksp_type = richardson
```

#### Quasistatic Poroelasticity

```{code-block} cfg
---
caption: PETSc options used for quasistatic poroelasticity in serial.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = ilu
```

```{code-block} cfg
---
caption: PETSc options used for quasistatic poroelasticity in parallel.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = gamg
mg_levels_pc_type = sor
mg_levels_ksp_type = richardson
```

### Monitoring

The monitoring options are enabled by default and provide a few lines of output per time step summarizing the operation of the linear and nonlinear solvers and time stepping.
Additional monitoring can be turned on using the user-specified options.

```{code-block} cfg
---
caption: Default PETSc options for monitoring.
---
# Turn off monitoring options (enabled by default).
[pylithapp.problem.petsc_defaults]
monitoring = False

# Corresponding PETSc options
[pylithapp.petsc]
ksp_converged_reason = true

snes_converged_reason = true
snes_monitor = true

ts_monitor = tru
ts_error_if_step_fails = true
```

### Testing

The options in the testing category are intended for use in internal testing.
These options help identify memory leaks in PETSc data structures and inconsistent back traces.

```{code-block} cfg
---
caption: Default PETSc options for testing.
---
# Turn on testing options (turned off by default).
[pylithapp.problem.petsc_defaults]
testing = True

# Corresponding PETSc options
[pylithapp.petsc]
malloc_dump = true
```

## User-Specified PETSc Options

{numref}`tab-petsc-options-monitor` shows the main monitoring options offered by PETSc.
When optimizing and troubleshooting solver options, we usually turn on all the monitoring.

```{table} Description of PETSc monitoring options
:name: tab-petsc-options-monitor
| Option | Description |
| :------| :-----------|
| `log_view` | Show logging objects and events. |
| `ts_monitor` | Show time-stepping progress. |
| `ksp_monitor` | Show preconditioned residual norm. |
| `ksp_view` | Show linear solver parameters. |
| `ksp_error_if_not_converged` | Generate an error if linear solver does not converge. |
| `ksp_converged_reason` | Indicate why iterating stopped in linear solve. |
| `snes_monitor` | Show residual norm for each nonlinear solve iteration. |
| `snes_view` | Show nonlinear solver parameters. |
| `snes_error_if_not_converged` | Generate an error if nonlinear solver does not converge. |
| `snes_converged_reason` | Indicate why iterating stopped in nonlinear solve. |
| `snes_linesearch_monitor` | Show line search information in nonlinear solve. |
```

### Solver Options

For most problems we use the GMRES method from Saad and Schultz for the linear solver; this is the linear solver PETSc uses as the default.
See [PETSc linear solver table](https://petsc.org/release/docs/manual/ksp/#tab-kspdefaults) for a list of PETSc options for linear solvers and preconditioners.

:::{tip}
It is important to keep in mind the resolution of the model and observations when setting solver tolerances.
For example, matching observations with an accuracy of 1.0 mm does not require solving the equations to an accuracy of 0.0001 mm.
:::

```{table} Summary of PETSc solver tolerances.
:name: tab-petsc-options-solver
| Option   | Description            |
| :--------- | :--------------------- |
| `ksp_rtol` | Stop iterating when the preconditioned KSP residual norm has this amount relative to its starting value. |
| `ksp_atol` | Stop iterating when the preconditioned KSP residual normal is smaller than this value. |
| `snes_rtol` | Stop iterating when the SNES residual norm has this amount relative to its starting value. |
| `snes_atol` | Stop iterating when the SNES residual normal is smaller than this value. |
```

% End of file
