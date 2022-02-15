(sec-user-petsc-settings)=
# PETSc Settings

PyLith relies on PETSc for the finite-element data structures, linear and nonlinear solvers, and time-stepping algorithms.
PETSc has its own object-oriented interface for specifying runtime options.
Instead of trying to maintain a Pyre interface to all of the PETSc options, we use a single `petsc` facility to collect all of the PETSc options and pass them to PETSc.

PETSc time-stepping options are discussed in {ref}`sec-user-problems-time-dependent`.

## Monitor and Logging Settings

{numref}`tab-petsc-options-monitor` shows the main monitoring options offered by PETSc.
When optimizing and troubleshooting solver settings, we usually turn on all the monitoring.
Our recommended settings for all simulations include:

```{code-block} cfg
---
caption: Recommended PETSc monitoring settings as set in a `cfg` file.
---
[pylithapp.petsc]
# Trigger errors if linear or nonlinear solver fails to converge.
ksp_error_if_not_converged = true
snes_error_if_not_converged = true

# Monitor converged reasons
ksp_converged_reason = true
snes_converged_reason = true

# Monitor time-stepping and nonlinear solver
ts_monitor = true
snes_monitor = true
snes_linesearch_monitor = true
```

```{table} Description of PETSc monitoring settings
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

## Solver Settings

For most problems we use the GMRES method from Saad and Schultz for the linear solver with solver tolerances around 1.0e-10.
When running large problems, we often raise the solver tolerances by one or two orders of magnitude to reduce runtime while still achieving suitable accuracy.

See [PETSc linear solver table](https://petsc.org/release/docs/manual/ksp/#tab-kspdefaults) for a list of PETSc options for linear solvers and preconditioners.

:::{tip}
It is important to keep in mind the resolution of the model and observations when setting solver tolerances.
For example, matching observations with an accuracy of 1.0 mm does not require solving the equations to an accuracy of 0.0001 mm.
:::

```{table} Recommended starting point for PETSc solver tolerances.
:name: tab-petsc-options-solver
| Property   | Value | Description            |
| :--------- | :--------: | :--------------------- |
| `ksp_rtol` | 1.0e-10 &nbsp;| Stop iterating when the preconditioned KSP residual norm has this amount relative to its starting value. |
| `ksp_atol` | 1.0e-12 | Stop iterating when the preconditioned KSP residual normal is smaller than this value. |
| `snes_rtol` | 1.0e-10 | Stop iterating when the SNES residual norm has this amount relative to its starting value. |
| `snes_atol` | 1.0e-10 | Stop iterating when the SNES residual normal is smaller than this value. |
```

### Small elasticity problems

When running small test problems (about 1k or less unknowns) it is very handy to use a robust preconditioner so that issues related to the boundary conditions, material parameters, etc. are more obvious.
We recommend using Incomplete (ILU) factorization.

```{code-block} cfg
---
caption: Recommended PETSc solver settings for small problems.
---
[pylithapp.petsc]
pc_type = ilu
ksp_type = gmres
```

### Medium elasticity problems

When running slightly larger problems (about 10k or less unknowns), the Additive Schwarz Method (ASM) using Incomplete LU (ILU) factorization preconditioner is usually more efficient.

```{code-block} cfg
---
caption: Recommended PETSc solver settings for medium problems
---
[pylithapp.petsc]
pc_type = asm
ksp_type = gmres
```

### Large elasticity problems without a fault}

Algebraic multigrid preconditioner usually works very well on elasticity problems.

```{code-block} cfg
---
caption: Recommended PETSc solver settings for solving elasticity problems without a fault
---
[pylithapp.petsc]
pc_type = ml
ksp_type = gmres
```

:::{important}
The ML algebraic multigrid preconditioner is only available if you build PETSc with the ML package. These features are included in the PyLith binary packages.
:::

### Large elasticity problems with a fault

The Lagrange multiplier solution subfield introduces a saddle point in the system of equations, so we use a Schur complement approach. These settings are available in `$PYLITH_DIR/share/settings/solver_fault_exact.cfg`.

```{code-block} cfg
---
caption: Recommended PETSc solver settings for solving elasticity problems with a fault.
---
[pylithapp.petsc]
pc_type = fieldsplit
pc_use_amat = true
pc_fieldsplit_type = schur
pc_fieldsplit_schur_factorization_type = full
pc_fieldsplit_dm_splits = true
fieldsplit_displacement_ksp_type = preonly
fieldsplit_displacement_pc_type = lu
fieldsplit_lagrange_multiplier_fault_pc_type = jacobi
fieldsplit_lagrange_multiplier_fault_ksp_type = gmres
fieldsplit_lagrange_multiplier_fault_ksp_rtol = 1.0e-11
fieldsplit_lagrange_multiplier_fault_ksp_converged_reason = true
```

:::{warning}
The split fields and algebraic multigrid preconditioning currently fails in problems with a nonzero null space.
This most often occurs when a problem contains multiple faults that extend through the entire domain and create subdomains without any Dirichlet boundary conditions.
The current workaround is to use the Additive Schwarz preconditioner without split fields.
See Section \vref{sec:Troubleshooting} for the error message encountered in this situation.
:::

### Incompressibility elasticity problems

The pressure solution subfield introduces a saddle point in the system of equations, so we again use a Schur complement approach.
This time we can use algebraic multigrid preconditioning on each block.
These settings are available in `$PYLITH_DIR/share/settings/solver_incompressible_elasticity.cfg`.

```{code-block} cfg
---
caption: Recommended PETSc solver settings for solving incompressible elasticity problems without a fault
---
[pylithapp.petsc]
pc_type = fieldsplit
pc_fieldsplit_type = schur
pc_fieldsplit_schur_fact_type = full
pc_fieldsplit_schur_precondition = full
fieldsplit_displacement_pc_type = lu
fieldsplit_pressure_pc_type = lu
```

% End of file
