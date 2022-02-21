# Governing Equations

This chapter presents the solution schemes we use for solving variations of the elasticity equation using the finite-element method.
In all of our derivations, we use the notation described in {ref}`tab:notation`.

```{table} Mathematical notation
:name: tab:notation
|    **Symbol**     | **Description**             |
| :---------------: | :-------------------------- |
|     $\vec{a}$     | Vector field a              |
|   $\mathbf{a}$    | Second order tensor field a |
|     $\vec{u}$     | Displacement vector field   |
|     $\vec{v}$     | Velocity vector field       |
|    $\vec{{d}}$    | Fault slip vector field     |
|     $\vec{f}$     | Body force vector field     |
|   $\vec{\tau}$    | Traction vector field       |
| $\mathbf{\sigma}$ | Stress tensor field         |
|     $\vec{n}$     | Normal vector field         |
|      $\rho$       | Mass density scalar field   |
```

:::{toctree}
elasticity-derivation.md
petsc-formulation.md
elasticity-infstrain/index.md
elasticity-infstrain-prescribedslip/index.md
incompressible-elasticity.md
poroelasticity/index.md
:::
