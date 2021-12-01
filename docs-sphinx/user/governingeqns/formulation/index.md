# Multiphysics Finite-Element Formulation

Within the PETSc solver framework, we want to solve a system of partial differential equations in which the weak form can be expressed as $F(t,s,\dot{s}) = G(t,s)$, $s(t_0) = s_0$, where $F$ and $G$ are vector functions, $t$ is time, and $s$ is the solution vector.

Using the finite-element method we manipulate the weak form of the system of equations involving a vector field $\vec{u}$ into integrals over the domain $\Omega$ matching the form,

$$
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{f}_0(t,s,\dot{s}) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : \boldsymbol{f}_1(t,s,\dot{s}) \, d\Omega =   \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{g}_0(t,s) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : \boldsymbol{g}_1(t,s) \, d\Omega,
$$ (eqn:problem:form)

where ${\vec{\psi}_\mathit{trial}^{u}}$ is the trial function for field $\vec{u}$, $\vec{f}_0$ and $\vec{g}_0$ are vectors, and $\boldsymbol{f}_1$ and $\boldsymbol{g}_1$ are tensors.
With multiple partial differential equations we will have multiple equations of this form, and the solution vector $s$, which we usually write as $\vec{s}$, will be composed of several different fields, such as displacement $\vec{u}$, velocity $\vec{v}$, pressure $p$, and temperature $T$.
Boundary conditions will also contribute similar terms with integrals over the corresponding boundaries.

For consistency with the PETSc time stepping formulation, we call $G(t,s)$ the RHS function and call $F(t,s,\dot{s})$ the LHS (or I) function.
Likewise, the Jacobian of $G(t,s)$ is the RHS Jacobian and the Jacobian of $F(t,s,\dot{s})$ is the LHS Jacobian.
Using a finite-element discretization we break up the domain and boundary integrals into sums over the cells and boundary faces/edges, respectively.
Using numerical quadrature those sums in turn involve sums over the values at the quadrature points with appropriate weights.
Thus, computation of the RHS function boils down to pointwise evaluation of $\vec{g}_0(t,s)$ and $\boldsymbol{g}_1(t,s)$, and computation of the LHS function boils down to pointwise evaluation of $\vec{f}_0(t,s,\dot{s})$ and $\boldsymbol{f}_1(t,s,\dot{s})$.

:::{toctree}
jacobian.md
petsc-notes.md
:::
