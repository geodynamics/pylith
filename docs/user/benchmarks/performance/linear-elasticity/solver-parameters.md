# Solver Parameters

We consider the solver parameters used in PyLith versions 3.x (Schur complement) and 4.x (variable point block Jacobi), including adjustments that improve performance.

## Schur Complement

In PyLith versions 3.0--4.0, the default solver parameters for linear elasticity with prescribed slip use a Schur complement approach ({numref}`listing:solver:parameters:fieldsplit:selfp`).
We use lower factorization (`pc_fieldsplit_schur_factorization_type = lower`) with the preconditioning for the Schur complement generated from an explicitly-assembled approximation (`pc_fieldsplit_schur_precondition = selfp`), $S_p=A_{11} − A_{10} (\mathit{diag}(A_{00}))^{-1} A_{01}$ \cite{petsc}.
Within the blocks, we use algebraic multigrid preconditioning with only one application of the preconditioner (`preonly`).

```{code-block} cfg
---
caption: Schur complement solver parameters (`fieldsplit (selfp)`).
name: listing:solver:parameters:fieldsplit:selfp
---
[pylithapp.problem.scales]
length_scale = 1.0*km
displacement_scale = 1.0*km

[pylithapp.problem.petsc_defaults]
solver = False

[pylithapp.petsc]
pc_type = fieldsplit
pc_use_amat = true
pc_fieldsplit_type = schur

pc_fieldsplit_schur_factorization_type = lower
pc_fieldsplit_schur_precondition = selfp
pc_fieldsplit_schur_scale = 1.0

fieldsplit_displacement_ksp_type = preonly
fieldsplit_displacement_pc_type = ml

fieldsplit_lagrange_multiplier_fault_ksp_type = preonly
fieldsplit_lagrange_multiplier_fault_pc_type = ml
```

## Variable Point Block Jacobi

In PyLith versions 4.1--4.2, the default solver parameters for linear elasticity with prescribed slip use a variable point block Jacobian preconditioner ({numref}`listing:solver:parameters:vpbjacobi:defaults`).
We reorder the degrees of freedom in the solution vector (`dm_reorder_section_type = cohesive`) to form small blocks with the degrees of freedom associated with the Lagrange multipliers at a point and the corresponding displacements.
We use the PETSc geometric algebraic multigrid preconditioer with variable point block Jacobi preconditioning at the fine scale (`mg_fine_pc_type = vpbjacobi`) and the  on the rest of the system of equations.
We use the default parameters for the geometric algebraic multigrid preconditioner.

```{code-block} cfg
---
caption: Variable point block Jacobi solver parameters (`vpbjacobi (defaults)`). The highlighted lines indicate the differences relative to the field split solver parameters.
name: listing:solver:parameters:vpbjacobi:defaults
emphasize-lines: 8-15
---
[pylithapp.problem.scales]
length_scale = 1.0*km
displacement_scale = 1.0*km

[pylithapp.problem.petsc_defaults]
solver = False

[pylithapp.petsc]
pc_type = gamg

dm_reorder_section = True
dm_reorder_section_type = cohesive

mg_fine_pc_type = vpbjacobi
ksp_gmres_restart = 100
```

### Tuning Fine-Scale Parameters

In PyLith version 5.0, we adjust the fine-scale preconditioner settings to improve solver performance (`mg_fine_ksp_max_it = 5`).

```{code-block} cfg
---
caption: Variable point block Jacobi solver parameters with adjustments at the fine scale of the multigrid preconditioner (`vpbjacobi (tunefine)`).
emphasize-lines: 16
---
[pylithapp.problem.scales]
length_scale = 1.0*km
displacement_scale = 1.0*km

[pylithapp.problem.petsc_defaults]
solver = False

[pylithapp.petsc]
pc_type = gamg

dm_reorder_section = True
dm_reorder_section_type = cohesive
mg_fine_pc_type = vpbjacobi
ksp_gmres_restart = 100

mg_fine_ksp_max_it = 5
```

### Tuning Fine- and Coarse-Scale Parameters

We also considered tuning at the parameters at the coarse levels of the multigrid preconditioner ({numref}`listing:solver:parameters:vpbjacobi:tunecoarse`).
These changes resulted in slight increases in the number of iterations required for the linear solver and a small increase in the runtime.

```{code-block} cfg
---
caption: Variable point block Jacobi solver parameters with adjustments at the fine and coarse scales of the multigrid preconditioner (`vpbjacobi (tunecoarse)`). Adjusting the parameters and the coarse scale did not improve performance.
emphasize-lines: 18-20
---
[pylithapp.problem.scales]
length_scale = 1.0*km
displacement_scale = 1.0*km

[pylithapp.problem.petsc_defaults]
solver = False

[pylithapp.petsc]
pc_type = gamg

dm_reorder_section = True
dm_reorder_section_type = cohesive
mg_fine_pc_type = vpbjacobi
ksp_gmres_restart = 100

mg_fine_ksp_max_it = 5

pc_gamg_coarse_eq_limit = 200
mg_levels_pc_type = pbjacobi
pc_gamg_agg_nsmooths = 3
```

### Length and Displacement Scales

In PyLith version 5.0, we also added independent specification displacement scale, whereas in previous versions the displacement scale matched the length scale.
Setting the displacement scale to the displacement corresponding to the average slip ({numref}`listing:solver:parameters:vpbjacobi:dispscale`) normalizes the displacement and Lagrange multiplier values in the solution vector independent of the length scale.
Consequently, the solver tolerances become independent of the discretization size and we can specify meaningful defaults (`ksp_atol=1.0e-7, snes_atol=4.0e-7`).

```{code-block} cfg
---
caption: Variable point block Jacobi solver parameters with adjustments at the fine scale of the multigrid preconditioner and separation of length and displacement scales (`vpbjacobi (dispscale)`).
emphasize-lines: 3
---
[pylithapp.problem.scales]
length_scale = 1.0*km
displacement_scale = 1.0*m

[pylithapp.problem.petsc_defaults]
solver = False

[pylithapp.petsc]
pc_type = gamg

dm_reorder_section = True
dm_reorder_section_type = cohesive
mg_fine_pc_type = vpbjacobi
ksp_gmres_restart = 100

mg_fine_ksp_max_it = 5
```

## Solver Tolerances

With the variable point block Jacobi preconditioner at the fine scale, for the scales to match in the blocks that include displacement and Lagrange multiplier terms, the length scale ($x_0$) must match the nominal discretization size ($h_0$).
The displacement and Lagrange multiplier values in the solution vector are normalized by the length scale, so the normalization depends on the discretization size.
To achieve a given level of accuracy in the displacement field relative to the scale of deformation requires setting the KSP and SNES absolute tolerances relative to the discretization size.

For example, to achieve an accuracy in the displacement field of 1.0e-7 relative to the nominal displacement, $u_\mathit{nominal}$, we use

\begin{gather}
\text{ksp_atol} = 1.0 \times 10^{-7} \ u_\mathit{nominal} / h_o, \\
\text{snes_atol} = 4.0 \times 10^{-7} \ u_\mathit{nominal} / h_o.
\end{gather}
