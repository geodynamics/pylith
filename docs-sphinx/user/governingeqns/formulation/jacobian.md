# Jacobian

The LHS Jacobian {math}`J_F = \frac{\partial F}{\partial s} + s_\mathit{tshift} \frac{\partial F}{\partial \dot{s}}` and the RHS Jacobian {math}`J_G = \frac{\partial G}{\partial s}`, where {math}`s_\mathit{tshift}` arises from the temporal discretization.
We put the Jacobians for each equation into the form:

```{math} J_F = \int_\Omega \vec{\psi}_\mathit{trial} \cdot \mathbf{J}_{f0}(t,s,\dot{s}) \cdot \vec{\psi}_\mathit{basis} + \vec{\psi}_\mathit{trial} \cdot \mathbf{J}_{f1}(t,s,\dot{s}) : \nabla \vec{\psi}_\mathit{basis}  + \nabla \vec{\psi}_\mathit{trial} : \mathbf{J}_{f2}(t,s,\dot{s}) \cdot \vec{\psi}_\mathit{basis} + \nabla \vec{\psi}_\mathit{trial} : \mathbf{J}_{f3}(t,s,\dot{s}) : \nabla \vec{\psi}_\mathit{basis} \  d\Omega
---
label: eqn:jacobian:form
---
```
```{math} J_G = \quad \int_\Omega \vec{\psi}_\mathit{trial} \cdot \mathbf{J}_{g0}(t,s) \cdot \vec{\psi}_\mathit{basis} + \vec{\psi}_\mathit{trial} \cdot \mathbf{J}_{g1}(t,s) : \nabla \vec{\psi}_\mathit{basis} + \nabla \vec{\psi}_\mathit{trial} : \mathbf{J}_{g2}(t,s) \cdot \vec{\psi}_\mathit{basis} + \nabla \vec{\psi}_\mathit{trial} : \mathbf{J}_{g3}(t,s) : \nabla \vec{\psi}_\mathit{basis} \  d\Omega,
```
where {math}`\vec{\psi}_\mathit{basis}` is a basis function.
Expressed in index notation the Jacobian coupling solution field components {math}`s_{i}` and {math}`s_{j}` is
```{math} J^{s_{i}s_{j}} = \int_\Omega \psi_{\mathit{trial},i} J_0^{s_{i}s_{j}} \psi_{\mathit{basis},j} + \psi_{\mathit{trial},i} J_1^{s_{j}s_{j}l} \psi_{\mathit{basis},j,l} + \psi_{\mathit{trial},i,k} J_2^{s_{i}s_{j}k} \psi_{\mathit{basis},j} + \psi_{\mathit{trial},i,k} J_3^{s_{i}s_{j}kl} \psi_{\mathit{basis},j,l} \  d\Omega,
---
label: eqn:jacobian:index:form
---
```
It is clear that the tensors {math}`J_0`, {math}`J_1`, {math}`J_2`, and {math}`J_3` have various sizes: {math}`J_0(n_i,n_j)`, {math}`J_1(n_i,n_j,d)`, {math}`J_2(n_i,n_j,d)`, {math}`J_3(n_i,n_j,d,d)`, where {math}`n_i` is the number of components in solution field {math}`s_{i}`, {math}`n_j` is the number of components in solution field {math}`s_{j}`, and {math}`d` is the spatial dimension.
Alternatively, expressed in discrete form, the Jacobian for the coupling between solution fields {math}`s_{i}` and {math}`s_{j}` is
```{math} J^{s_{i}s_{j}} = J_{0}^{s_{i}s_{j}} + J_{1}^{s_{i}s_{j}} B + B^T J_{2}^{s_{i}s_{j}} + B^T J_{3}^{s_{i}s_{j}} B,
---
label: eqn:jacobian:discrete:form
```
where {math}`B` is a matrix of the derivatives of the basis functions and {math}`B^T` is a matrix of the derivatives of the trial functions.

:::{important}
See <https://www.mcs.anl.gov/petsc/petsc-master/docs/manualpages/FE/PetscFEIntegrateJacobian.html> for the ordering of indices in the Jacobian pointwise functions.
:::
