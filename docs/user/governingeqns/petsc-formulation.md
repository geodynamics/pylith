(sec-user-petsc-fe-formulation)=
# Finite-Element Formulation with PETSc

Within the PETSc solver framework, we want to solve a system of partial differential equations in which the weak form can be expressed as $F(t,s,\dot{s}) = G(t,s)$, $s(t_0) = s_0$, where $F$ and $G$ are vector functions, $t$ is time, and $s$ is the solution vector.

Using the finite-element method[^1] we manipulate the weak form of the system of equations involving a vector field $\vec{u}$ into integrals over the domain $\Omega$ matching the form,

```{math}
:label: eqn:problem:form
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{f}_0(t,s,\dot{s}) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : \boldsymbol{f}_1(t,s,\dot{s}) \, d\Omega =   \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{g}_0(t,s) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : \boldsymbol{g}_1(t,s) \, d\Omega,
```

where ${\vec{\psi}_\mathit{trial}^{u}}$ is the trial function for field $\vec{u}$, $\vec{f}_0$ and $\vec{g}_0$ are vectors, and $\boldsymbol{f}_1$ and $\boldsymbol{g}_1$ are tensors.
With multiple partial differential equations we will have multiple equations of this form, and the solution vector $s$, which we usually write as $\vec{s}$, will be composed of several different fields, such as displacement $\vec{u}$, velocity $\vec{v}$, pressure $p$, and temperature $T$.
Boundary conditions will also contribute similar terms with integrals over the corresponding boundaries.

For consistency with the PETSc time stepping formulation, we call $G(t,s)$ the RHS function and call $F(t,s,\dot{s})$ the LHS (or I) function.
Likewise, the Jacobian of $G(t,s)$ is the RHS Jacobian and the Jacobian of $F(t,s,\dot{s})$ is the LHS Jacobian.
Using a finite-element discretization we break up the domain and boundary integrals into sums over the cells and boundary faces/edges, respectively.
Using numerical quadrature those sums in turn involve sums over the values at the quadrature points with appropriate weights.
Thus, computation of the RHS function boils down to pointwise evaluation of $\vec{g}_0(t,s)$ and $\boldsymbol{g}_1(t,s)$, and computation of the LHS function boils down to pointwise evaluation of $\vec{f}_0(t,s,\dot{s})$ and $\boldsymbol{f}_1(t,s,\dot{s})$.

## Jacobian

The LHS Jacobian $J_F = \frac{\partial F}{\partial s} + s_\mathit{tshift} \frac{\partial F}{\partial \dot{s}}$ and the RHS Jacobian $J_G = \frac{\partial G}{\partial s}$, where $s_\mathit{tshift}$ arises from the temporal discretization. We put the Jacobians for each equation into the form:

```{math}
:label: eqn:jacobian:form
\begin{aligned}
  J_F &= \int_\Omega {\vec{\psi}_\mathit{trial}^{}}\cdot \boldsymbol{J}_{f0}(t,s,\dot{s}) \cdot {\vec{\psi}_\mathit{basis}^{}} + {\vec{\psi}_\mathit{trial}^{}}\cdot \boldsymbol{J}_{f1}(t,s,\dot{s}) : \nabla {\vec{\psi}_\mathit{basis}^{}} + \nabla {\vec{\psi}_\mathit{trial}^{}}: \boldsymbol{J}_{f2}(t,s,\dot{s}) \cdot {\vec{\psi}_\mathit{basis}^{}} + \nabla {\vec{\psi}_\mathit{trial}^{}}: \boldsymbol{J}_{f3}(t,s,\dot{s}) : \nabla {\vec{\psi}_\mathit{basis}^{}}\, d\Omega \\
%
  J_G &= \int_\Omega {\vec{\psi}_\mathit{trial}^{}}\cdot \boldsymbol{J}_{g0}(t,s) \cdot {\vec{\psi}_\mathit{basis}^{}} + {\vec{\psi}_\mathit{trial}^{}}\cdot \boldsymbol{J}_{g1}(t,s) : \nabla {\vec{\psi}_\mathit{basis}^{}} + \nabla {\vec{\psi}_\mathit{trial}^{}}: \boldsymbol{J}_{g2}(t,s) \cdot {\vec{\psi}_\mathit{basis}^{}} + \nabla {\vec{\psi}_\mathit{trial}^{}}: \boldsymbol{J}_{g3}(t,s) : \nabla {\vec{\psi}_\mathit{basis}^{}}\, d\Omega,
\end{aligned}
```

where ${\vec{\psi}_\mathit{basis}^{}}$ is a basis function.
Expressed in index notation the Jacobian coupling solution field components $s_i$ and $s_j$ is

```{math}
:label: (eqn:jacobian:index:form)
J^{s_is_j} = \int_\Omega {\psi_\mathit{trial}^{}}_i J_0^{s_is_j} {\psi_\mathit{basis}^{}}_j + {\psi_\mathit{trial}^{}}_i
J_1^{s_js_jl}
{\psi_\mathit{basis}^{}}_{j,l} + {\psi_\mathit{trial}^{}}_{i,k} J_2^{s_is_jk} {\psi_\mathit{basis}^{}}_j + {\psi_\mathit{trial}^{}}_{i,k}
J_3^{s_is_jkl}
{\psi_\mathit{basis}^{}}_{j,l} \, d\Omega,
```

It is clear that the tensors $J_0$, $J_1$, $J_2$, and $J_3$ have various sizes: $J_0(n_i,n_j)$, $J_1(n_i,n_j,d)$, $J_2(n_i,n_j,d)$, $J_3(n_i,n_j,d,d)$, where $n_i$ is the number of components in solution field $s_i$, $n_j$ is the number of components in solution field $s_j$, and $d$ is the spatial dimension.
Alternatively, expressed in discrete form, the Jacobian for the coupling between solution fields $s_i$ and $s_j$ is

```{math}
:label: eqn:jacobian:discrete:form
  J^{s_is_j} = J_{0}^{s_is_j} + J_{1}^{s_is_j} B + B^T J_{2}^{s_is_j} + B^T J_{3}^{s_is_j} B,
```

where $B$ is a matrix of the derivatives of the basis functions and $B^T$ is a matrix of the derivatives of the trial functions.

:::{important}
See <https://www.mcs.anl.gov/petsc/petsc-master/docs/manualpages/FE/PetscFEIntegrateJacobian.html> for the ordering of indices in the Jacobian pointwise functions.
:::

## PETSc TS Notes

### Explicit Time Stepping

Explicit time stepping with the PETSc TS requires $F(t,s,\dot{s}) = \dot{s}$.
* We do not specify the functions $\vec{f}_0(t,s,\dot{s})$ and $\boldsymbol{f}_1(t,s,\dot{s})$ because the PETSc TS will assume $F(t,s,\dot{s}) = \dot{s}$ if no LHS (or I) function is given.
* The PETSc TS will verify that the LHS (or I) function is null.
* We also do not specify $J_F$ or $J_G$.
* This leaves us with only needing to specify $\vec{g}_0(t,s)$ and $\boldsymbol{g}_1(t,s)$.

For explicit time stepping with the PETSc TS, we need $F(t,s,\dot{s}) = \dot{s}$.
Using a finite-element formulation for elastodynamics, $F(t,s,\dot{s})$ generally involves integrals of the inertia over the domain.
It is tempting to simply move these terms to the RHS, but that introduces inertial terms into the boundary conditions, which makes them less intuitive.
Instead, we transform our equation into the form $\dot{s} = G^*(t,s)$ where $G^*(t,s) = M^{-1} G(t,s)$.
We take $M$ to be a lumped (diagonal) matrix, so that $M^{-1}$ is trivial to compute.
In computing the RHS function, $G^*(t,s)$, we compute $G(t,s)$, then compute $M$ and $M^{-1}$, and then $M^{-1}G(t,s)$.
For the elasticity equation with inertial terms, $M$ contains the mass matrix.

### Implicit Time Stepping

The LHS (or I) function is associated with implicit time-stepping.
When using implicit time-stepping, we place all of the terms on the LHS.
Even though placing all of the terms on the LHS sometimes requires different pointwise functions for implicit and explicit time stepping, it minimizes the number of pointwise functions needed for implicit time stepping.
If no RHS function is given, then the PETSc TS assumes $G(t,s) = 0$, so we only need to specify $F(t,s,\dot{s})$ and $J_f$.

### Implicit-Explicit Time Stepping

For implicit-explicit time stepping algorithms, the equations integrated with explicit time stepping have $\dot{s}$ as the LHS function, and the equations integrated with implicit time stepping have 0 as the RHS function.

[^1]: Some resources for learning about the finite-element method include {cite:t}`Zienkiewicz:Taylor:2000`, {cite:t}`Taylor:2003`,  {cite:t}`Farrell:FEM:2021`, and {cite:t}`Cotter:Ham:FEM:2023`.
