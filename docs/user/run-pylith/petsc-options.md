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

:::{note}
New in v3.0.0
:::

We separate the defaults into a few categories to make it easy to select desired options.

:monitors: Options for basic monitoring of the solver;
:initial_guess: Options for improving initial guesses at each time step;
:collective_io: Options for collective input/output; and
:testing: Options used in testing.
:solver: Options for the preconditioner and solver;
:parallel: Options used when running in parallel (can be used in serial as well);

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

### Monitoring Options

The monitoring options are enabled by default and provide a few lines of output per time step summarizing the operation of the linear and nonlinear solvers and time stepping.
Additional monitoring can be turned on using the user-specified options.

```{table} Description of PETSc monitoring options
:name: tab-petsc-options-monitor
| Option                        | Description                                              |
| :---------------------------- | :------------------------------------------------------- |
| `ts_monitor`                  | Show time-stepping progress.                             |
| `ksp_monitor`                 | Show preconditioned residual norm.                       |
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
ts_monitor = True
ksp_monitor = True
snes_monitor = True

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
Additionally, PyLith specifies general options related to the solver tolerances and triggering errors if the linear or nonlinear solver fails to converge.
The different sets of defaults are detailed in the following code blocks.

:::{warning}
When running in parallel in cases in which a fault face is split across processes, the current solver settings will result in a diverged solution.
We attempt to prevent this from happening by specifying a penalty for splitting across the fault; however, sometimes the partitioner will still split a fault face across processes.
In cases in which the default settings fail and the solver diverges, you can fall back to the previous settings by using the field split preconditioner in `share/settings/solver_elasticity_fault_fieldsplit.cfg`.
Simply add this `.cfg` file to your command line options.

This issue will go away once we implement parallel mesh loading.

:::{note}
If you do use the field split fall back, you need to be aware that it also has deficiencies.
The split fields and algebraic multigrid preconditioning currently fails in problems with a nonzero null space.
This most often occurs when a problem contains multiple faults that extend through the entire domain and create subdomains without any Dirichlet boundary conditions.
The workaround is to use the `ilu` preconditioner.
However, it only works in serial.
An alternative is to use the `asm` preconditioner (Additive Schwarz) which works in parallel and serial.
:::

:::

#### Solver tolerances

PyLith will set default solver tolerances wheneve the solver defaults are enabled.

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
ksp_rtol = 1.0e-12
ksp_atol = 1.0e-12

snes_rtol = 1.0e-12
snes_atol = 1.0e-9
```

#### Quasistatic Elasticity

```{code-block} cfg
---
caption: PETSc options used for quasistatic elasticity without a fault. These settings are also in `share/settings/solver_elasticity.cfg`.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = gamg
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

dm_reorder_section = True
dm_reorder_section_type = cohesive
mg_fine_pc_type = vpbjacobi
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
caption: PETSc options used for quasistatic incompressible elasticity in parallel. These settings are also in `share/settings/solver_incompressible_elasticity.cfg`.
---
[pylithapp.petsc]
ts_type = beuler
pc_type = fieldsplit
pc_fieldsplit_type = schur

pc_fieldsplit_schur_factorization_type = full
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
fieldsplit_displacement_mg_fine_pc_type = vpbjacobi

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

