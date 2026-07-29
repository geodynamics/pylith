(sec-user-run-pylith-petsc-options)=
# PETSc Options

PyLith relies on PETSc for finite-element data structures, linear and nonlinear solvers, and time-stepping algorithms.
PETSc has its own object-oriented interface for specifying runtime options.
PyLith provides two mechanisms for passing options to PETSc:

1. Default values that are controlled by a few flags in the PyLith `Problem` `petsc_defaults` facility, and
2. PyLith `petsc` facility which passes options directly to PETSc.

When using the default values, PyLith selects solver and preconditioner options based on the governing equations.
In most cases these default values will give good performance and you do not need to specify any PETSc options.
User-specified values always take precedence over the default values.

## Default PETSc Options

*New in v3.0.0.*

We separate the defaults into a few categories to make it easy to select desired options.

:monitors: Options for basic monitoring of the solver;
:initial_guess: Options for improving initial guesses at each time step;
:adaptive_time_stepping: Options for adaptive time stepping
:collective_io: Options for collective input/output; and
:testing: Options used in testing.
:solver: Options for the preconditioner and solver;
:parallel: Options used when running in parallel (can be used in serial as well);

:::{seealso}
See [`PetscDefaults` Component](../components/utils/PetscDefaults.md) for more information about the Pyre interface for specifying default PETSc options.
:::

### Monitoring Options

The monitoring options are enabled by default and provide a few lines of output per time step summarizing the operation of the linear and nonlinear solvers and time stepping.
Additional monitoring can be turned on using the user-specified options.

```{table} Description of PETSc monitoring options
:name: tab-petsc-options-monitor
| Option                        | Description                                              |
| :---------------------------- | :------------------------------------------------------- |
| `ts_monitor`                  | Show time-stepping progress.                             |
| `ksp_monitor`                 | Show preconditioned residual norm.                       |
| `ksp_monitor_true_residual`   | Show preconditioned and true residual norm.              |
| `ksp_error_if_not_converged`  | Generate an error if linear solver does not converge.    |
| `ksp_converged_reason`        | Indicate why iterating stopped in linear solve.          |
| `snes_monitor`                | Show residual norm for each nonlinear solve iteration.   |
| `snes_error_if_not_converged` | Generate an error if nonlinear solver does not converge. |
| `snes_converged_reason`       | Indicate why iterating stopped in nonlinear solve.       |
| `snes_linesearch_monitor`     | Show line search information in nonlinear solve.         |
```

```{code-block} cfg
---
caption: Default PETSc options for monitoring the solver. These settings are also in `share/settings/petsc_monitoring.cfg`.
---
[pylithapp.petsc]
# Turn on TS, KSP, and SNES monitors
ts_monitor =
ksp_monitor =
snes_monitor =

ksp_converged_reason = True
snes_converged_reason = True

# Trigger error if linear or nonlinear solvers fail to converge
ts_error_if_step_fails = True
ksp_error_if_not_converged = True
snes_error_if_not_converged = True
```

### Initial Guess Options

Improve initial guesses of the solution when time stepping.

```{code-block} cfg
---
caption: Default PETSc options for initial guesses of the solution. These settings are also in `share/settings/petsc_initial_guess.cfg`.
---
[pylithapp.petsc]
ksp_guess_type = pod
ksp_guess_pod_size = 8
```

### Adaptive Time-Stepping Options

*New in v5.0.0.*

Adaptive time stepping (not enabled by default) adjusts the time step based upon the rate of deformation.
We use the [`basic` PETSc adaptive time-stepping algorithm](https://petsc.org/release/manual/ts/#tab-adaptors).
Tolerances (`ts_atol` and `ts_rtol`) of around 0.05 work reasonably well in many simulations with viscoelasticity or poroelasticity.
Increasing the tolerances leds to relatively larger time steps and smaller tolerances lead to small time steps.
The `ts_adapt_clip`, `ts_adapt_safety` and `ts_adapt_reject_safety` parameters also influence the time step selection.
The [PETSc TS documentation](https://petsc.org/release/manual/ts/#error-control-via-variable-time-stepping) provides more details about how time steps are selected.

```{code-block} cfg
---
caption: Default PETSc options for adaptive time stepping. These settings are also in `share/settings/petsc_tsadaptive.cfg]`.
---
[pylithapp.petsc]
ts_adapt_type = basic
ts_adapt_safety = 0.2
ts_adapt_reject_safety = 0.1
ts_atol = 0.05
ts_rtol = 0.05

ts_adapt_monitor = true
```

### Collective I/O Options

The collective input and output options are enabled by default and turn on HDF5 collective output.
We use parallel HDF5 implementation, which in turn relies on MPI IO.
Many MPI IO implementations require collective input and output to be enabled for parallel HDF5 output even if only one process is being used.

```{code-block} cfg
---
caption: Default PETSc options for collective output. These settings are also in `share/settings/petsc_collectiveio.cfg`.
---
[pylithapp.petsc]
viewer_hdf5_collective = True
```

### Testing Options

The options in the testing category are intended for use in internal testing.
These options help identify memory leaks in PETSc data structures and inconsistent back traces.

```{code-block} cfg
---
caption: Default PETSc options for testing. These settings are also in `share/settings/petsc_testing.cfg`.
---
[pylithapp.petsc]
malloc_dump = True
```

### Solver Options

The solver options are enabled by default.
PyLith selects options based on the governing equation, formulation, presence of a fault, and whether the simulation is running in parallel.
In some cases the solver settings for running in parallel are different than those for running in serial; in such cases, the settings for running in parallel often given give comparable or better performance.
If you have a moderate or large simulation, you should enable the parallel settings.
Additionally, PyLith specifies default solver tolerances and options for triggering errors if the linear or nonlinear solver fails to converge.

#### Solver tolerances

PyLith will set default solver tolerances whenever the solver defaults are enabled.

We set the KSP and SNES relative tolerances (`ksp_rtol` and `snes_rtol`) to extremely small values (1.0e-14) so that convergence is never triggered by the relative criterion; these tolerances specify the residual reduction relative to the initial residual, and setting them this small effectively disables them in favor of the absolute tolerances described below.

The SNES absolute tolerance (`snes_atol`) is the true residual norm we consider sufficient for convergence. We use a default value of 5.0e-7, which yields displacements and stresses accurate to about that fraction of their respective characteristic scales.
This provides ample accuracy for most crustal deformation problems and is also small enough to be meaningful in our full-scale tests.

The KSP absolute tolerance (`ksp_atol`) controls when the linear solver declares convergence, but it is evaluated on the *preconditioned* residual, not the true residual, so the two do not match directly.
Under good conditions — an appropriate nondimensionalization, an effective preconditioner, and a mesh free of badly distorted cells — the preconditioned residual is typically within about an order of magnitude of the true residual.
In typical real-world simulations, this factor can grow to about 100.

To account for these variable conditions, we set `ksp_atol` to be smaller than `snes_atol` by a factor of 5, ensuring the linear solver converges tightly enough to satisfy the SNES convergence criterion without performing excessive iterations.
If preconditioned and true residual norms diverge by more than expected — e.g., due to a poorly scaled problem or a badly distorted mesh — you should try to improve the mesh quality and then assess convergence by inspecting the true residual (`ksp_monitor_true_residual`) and adjust the `ksp_atol` if necessary.

```{table} Summary of PETSc solver tolerances.
:name: tab-petsc-options-solver
| Option      | Description                                                                                              |
| :---------- | :------------------------------------------------------------------------------------------------------- |
| `ksp_rtol`  | Stop iterating when the preconditioned KSP residual norm has this amount relative to its starting value. |
| `ksp_atol`  | Stop iterating when the preconditioned KSP residual normal is smaller than this value.                   |
| `snes_rtol` | Stop iterating when the SNES residual norm has this amount relative to its starting value.               |
| `snes_atol` | Stop iterating when the SNES residual normal is smaller than this value.                                 |
```

```{code-block} cfg
---
caption: PETSc solver tolerances used whenever the solver defaults are enabled. These settings are also in `share/settings/petsc_solver_tolerances.cfg`.
---
[pylithapp.petsc]
ksp_rtol = 1.0e-14
ksp_atol = 1.0e-7

snes_rtol = 1.0e-14
snes_atol = 5.0e-7
```

#### Quasistatic Elasticity

```{code-block} cfg
---
caption: PETSc options used for quasistatic elasticity without a fault. These settings are also in `share/settings/solver_elasticity.cfg`.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = gamg
mg_fine_ksp_max_it = 5
```

The Lagrange multiplier corresponding to the tractions on the fault introduces a saddle point in the system of equations.
We could use a Schur complement approach, but we have found that grouping the degrees of freedom on each side of the fault into blocks and using a variable point-block Jacobi preconditioner provides better results; the number of iterations remains nearly constant with increased problem size and the overall solution time is low.

```{code-block} cfg
---
caption: PETSc options used for quasistatic elasticity with a fault. These settings are also in `share/settings/solver_elasticity_fault.cfg`.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = gamg
ksp_gmres_restart = 100

dm_reorder_section = True
dm_reorder_section_type = cohesive
mg_fine_pc_type = vpbjacobi
mg_fine_ksp_max_it = 5
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

pc_fieldsplit_schur_factorization_type = lower
pc_fieldsplit_schur_precondition = full

fieldsplit_displacement_pc_type = lu
fieldsplit_pressure_pc_type = lu
```

```{code-block} cfg
---
caption: PETSc options used for quasistatic incompressible elasticity in parallel. These settings are also in `share/settings/solver_incompressible_elasticity.cfg`.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_fieldsplit_type = schur

pc_fieldsplit_schur_factorization_type = lower
pc_fieldsplit_schur_precondition = full

fieldsplit_displacement_pc_type = ml
fieldsplit_pressure_pc_type = bjacobi
```

#### Quasistatic Poroelasticity

```{code-block} cfg
---
caption: PETSc options used for quasistatic poroelasticity. These settings are also in `share/settings/solver_poroelasticity.cfg`.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_fieldsplit_type = multiplicative
pc_fieldsplit_0_fields = 2
pc_fieldsplit_1_fields = 1
pc_fieldsplit_2_fields = 0

fieldsplit_trace_strain_pc_type = bjacobi
fieldsplit_pressure_pc_type = bjacobi

fieldsplit_displacement_pc_type = ml
fieldsplit_displacement_ksp_type = gmres
```

```{code-block} cfg
---
caption: PETSc options used for quasistatic poroelasticity with a porosity state variable. These settings are also in `share/settings/solver_poroelasticity_statevars.cfg`.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_fieldsplit_type = multiplicative
pc_fieldsplit_0_fields = 2
pc_fieldsplit_1_fields = 1
pc_fieldsplit_2_fields = 0
pc_fieldsplit_3_fields = 3
pc_fieldsplit_4_fields = 4
pc_fieldsplit_5_fields = 5

fieldsplit_trace_strain_pc_type = bjacobi
fieldsplit_pressure_pc_type = bjacobi

fieldsplit_displacement_pc_type = ml
fieldsplit_displacement_ksp_type = gmres

fieldsplit_velocity_pc_type = bjacobi
fieldsplit_pressure_t_pc_type = bjacobi
fieldsplit_trace_strain_t_pc_type = bjacobi
```

```{code-block} cfg
---
caption: PETSc options used for quasistatic poroelasticity with a fault. These settings are also in `share/settings/solver_poroelasticity_fault.cfg`.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_fieldsplit_type = multiplicative
pc_fieldsplit_0_fields = 2
pc_fieldsplit_1_fields = 1
pc_fieldsplit_2_fields = 0
pc_fieldsplit_3_fields = 3

fieldsplit_trace_strain_pc_type = bjacobi
fieldsplit_pressure_pc_type = bjacobi

fieldsplit_displacement_pc_type = ml
fieldsplit_displacement_ksp_type = gmres
fieldsplit_displacement_ksp_gmres_restart = 100
fieldsplit_displacement_mg_fine_pc_type = vpbjacobi
fieldsplit_displacement_mg_fine_ksp_max_it = 5

fieldsplit_lagrange_multiplier_fault_pc_type = jacobi

dm_reorder_section = True
dm_reorder_section_type = cohesive
```

## User-Specified PETSc Options

PETSc provides a few options for understanding solver settings and performance.
The options `ksp_view` and `snes_view` control the display of the details of the linear and nonlinear solvers, respectively.
The `log_view` option controls output of logging information at the end of a simulation.

```{table} Description of PETSc logging and solving viewing options
:name: tab-petsc-options-view
| Option      | Description                       |
| :---------- | :-------------------------------- |
| `log_view`  | Show logging objects and events.  |
| `ksp_view`  | Show linear solver parameters.    |
| `snes_view` | Show nonlinear solver parameters. |
```

```{code-block} cfg
---
caption: Examples of using the logging and viewing options.
---
[pylithapp.petsc]

# Show detailed information about the linear solver.
ksp_view = True

# Show detailed information about the nonlinear solver.
snes_view = True

# Write logging information to stdout.
log_view = True

# Write logging information to an ASCII file `pylith_log.txt`.
log_view = ascii:pylith_log.txt
```

### Customizing Solver Options

For most problems we use the GMRES method from {cite:t}`Saad:Schultz:1986` for the linear solver; this is the linear solver PETSc uses as the default.
See [PETSc linear solver table](https://petsc.org/release/docs/manual/ksp/#tab-kspdefaults) for a list of PETSc options for linear solvers and preconditioners.

:::{tip}
It is important to keep in mind the resolution of the model and observations when setting solver tolerances.
For example, matching observations with an accuracy of 1.0 mm does not require solving the equations to an accuracy of 0.0001 mm.
:::

