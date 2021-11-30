# Multiphysics Finite-Element Formulation

Within the PETSc solver framework, we want to solve a system of partial differential equations in which the weak form can be expressed as {math}`F(t,s,\dot{s}) = G(t,s)`, {math}`s(t_0) = s_0`, where {math}`F` and {math}`G` are vector functions, {math}`t` is time, and {math}`s` is the solution vector.

Using the finite-element method we manipulate the weak form of the system of equations involving a vector field {math}`\vec{u}` into integrals over the domain {math}`\Omega` matching the form,
```{math} \label{eqn:problem:form} \int_\Omega \vec{\psi}_\mathit{trial}^{u} \cdot \vec{f}_0(t,s,\dot{s}) + \nabla \vec{\psi}_\mathit{trial}^{u} : \mathbf{f} _1(t,s,\dot{s}) \ d\Omega =   \int_\Omega \vec{\psi}_\mathit{trial}^{u} \cdot \vec{g}_0(t,s) + \nabla \vec{\psi}_\mathit{trial}^{u} : \mathbf{g}_1(t,s) \ d\Omega,
---
label: eqn:problem:form
---
```
where {math}`\vec{\psi}_\mathit{trial}^{u}` is the trial function for field {math}`\vec{u}`, {math}`\vec{f}_0` and {math}`\vec{g}_0` are vectors, and {math}`\mathbf{f}_1` and {math}`\mathbf{g}_1` are tensors.
With multiple partial differential equations we will have multiple equations of this form, and the solution vector {math}`s`, which we usually write as {math}`\vec{s}`, will be composed of several different fields, such as displacement {math}`\vec{u}`, velocity {math}`\vec{v}`, pressure {math}`p`, and temperature {math}`T`.
Boundary conditions will also contribute similar terms with integrals over the corresponding boundaries.

For consistency with the PETSc time stepping formulation, we call {math}`G(t,s)` the RHS function and call {math}`F(t,s,\dot{s})` the LHS (or I) function.
Likewise, the Jacobian of {math}`G(t,s)` is the RHS Jacobian and the Jacobian of {math}`F(t,s,\dot{s})` is the LHS Jacobian.
Using a finite-element discretization we break up the domain and boundary integrals into sums over the cells and boundary faces/edges, respectively.
Using numerical quadrature those sums in turn involve sums over the values at the quadrature points with appropriate weights.
Thus, computation of the RHS function boils down to pointwise evaluation of {math}`\vec{g}_0(t,s)` and {math}`\mathbf{g}_1(t,s)`, and computation of the LHS function boils down to pointwise evaluation of {math}`\vec{f}_0(t,s,\dot{s})` and {math}`\mathbf{f}_1(t,s,\dot{s})`.

:::{toctree}
jacobian.md
petsc-notes.md
:::
