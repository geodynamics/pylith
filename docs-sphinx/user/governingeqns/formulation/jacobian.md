# Jacobian

The LHS Jacobian $J_F = \frac{\partial F}{\partial s} + s_\mathit{tshift} \frac{\partial F}{\partial \dot{s}}$ and the RHS Jacobian $J_G = \frac{\partial G}{\partial s}$, where $s_\mathit{tshift}$ arises from the temporal discretization. We put the Jacobians for each equation into the form:

$$
\begin{aligned}
  J_F &= \int_\Omega {\vec{\psi}_\mathit{trial}^{}}\cdot \boldsymbol{J}_{f0}(t,s,\dot{s}) \cdot {\vec{\psi}_\mathit{basis}^{}} + {\vec{\psi}_\mathit{trial}^{}}\cdot \boldsymbol{J}_{f1}(t,s,\dot{s}) : \nabla {\vec{\psi}_\mathit{basis}^{}} + \nabla {\vec{\psi}_\mathit{trial}^{}}: \boldsymbol{J}_{f2}(t,s,\dot{s}) \cdot {\vec{\psi}_\mathit{basis}^{}} + \nabla {\vec{\psi}_\mathit{trial}^{}}: \boldsymbol{J}_{f3}(t,s,\dot{s}) : \nabla {\vec{\psi}_\mathit{basis}^{}}\, d\Omega \\
  J_G &= \quad \int_\Omega {\vec{\psi}_\mathit{trial}^{}}\cdot \boldsymbol{J}_{g0}(t,s) \cdot {\vec{\psi}_\mathit{basis}^{}}   + {\vec{\psi}_\mathit{trial}^{}}\cdot \boldsymbol{J}_{g1}(t,s) : \nabla {\vec{\psi}_\mathit{basis}^{}}   + \nabla {\vec{\psi}_\mathit{trial}^{}}: \boldsymbol{J}_{g2}(t,s) \cdot {\vec{\psi}_\mathit{basis}^{}} + \nabla {\vec{\psi}_\mathit{trial}^{}}: \boldsymbol{J}_{g3}(t,s) : \nabla {\vec{\psi}_\mathit{basis}^{}}\, d\Omega,
\end{aligned}
$$ (eqn:jacobian:form)

where ${\vec{\psi}_\mathit{basis}^{}}$ is a basis function.
Expressed in index notation the Jacobian coupling solution field components $s_i$ and $s_j$ is

$$
J^{s_is_j} = \int_\Omega {\psi_\mathit{trial}^{}}_i J_0^{s_is_j} {\psi_\mathit{basis}^{}}_j + {\psi_\mathit{trial}^{}}_i
J_1^{s_js_jl}
{\psi_\mathit{basis}^{}}_{j,l} + {\psi_\mathit{trial}^{}}_{i,k} J_2^{s_is_jk} {\psi_\mathit{basis}^{}}_j + {\psi_\mathit{trial}^{}}_{i,k}
J_3^{s_is_jkl}
{\psi_\mathit{basis}^{}}_{j,l} \, d\Omega,
$$ (eqn:jacobian:index:form)

It is clear that the tensors $J_0$, $J_1$, $J_2$, and $J_3$ have various sizes: $J_0(n_i,n_j)$, $J_1(n_i,n_j,d)$, $J_2(n_i,n_j,d)$, $J_3(n_i,n_j,d,d)$, where $n_i$ is the number of components in solution field $s_i$, $n_j$ is the number of components in solution field $s_j$, and $d$ is the spatial dimension.
Alternatively, expressed in discrete form, the Jacobian for the coupling between solution fields $s_i$ and $s_j$ is

$$
  J^{s_is_j} = J_{0}^{s_is_j} + J_{1}^{s_is_j} B + B^T J_{2}^{s_is_j} + B^T J_{3}^{s_is_j} B,
$$ (eqn:jacobian:discrete:form)

where $B$ is a matrix of the derivatives of the basis functions and $B^T$ is a matrix of the derivatives of the trial functions.

:::{important}
See <https://www.mcs.anl.gov/petsc/petsc-master/docs/manualpages/FE/PetscFEIntegrateJacobian.html> for the ordering of indices in the Jacobian pointwise functions.
:::
