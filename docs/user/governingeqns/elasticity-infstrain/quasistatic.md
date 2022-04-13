# Quastistatic

If we neglect the inertial term ($\rho \frac{\partial \vec{v}}{\partial t} \approx \vec{0}$), then time dependence only arises from history-dependent constitutive equations and boundary conditions.
Our solution vector is the displacement vector and the elasticity equation reduces to
%
```{math}
:label: eqn:elasticity:strong:form:quasistatic
\begin{gather}
\vec{f}(\vec{x},t) + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u}) = \vec{0} \text{ in }\Omega, \\
%
\boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on }\Gamma_\tau, \\
%
\vec{u} = \vec{u}_0(\vec{x},t) \text{ on }\Gamma_u.
\end{gather}
```
%
Because we will use implicit time stepping, we place all of the terms in the elasticity equation on the LHS.
We create the weak form by taking the dot product with the trial function ${\vec{\psi}_\mathit{trial}^{u}}$ and integrating over the domain:
%
\begin{equation}
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \left( \vec{f}(t) + \boldsymbol{\nabla}\cdot \boldsymbol{\sigma} (\vec{u}) \right) \, d\Omega = 0.
\end{equation}
%
Using the divergence theorem and incorporating the Neumann boundary conditions, we have
%
\begin{equation}
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\boldsymbol{\sigma}(\vec{u}) \, d\Omega  + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma = 0
\end{equation}
%
## Residual Pointwise Functions

Identifying $F(t,s,\dot{s})$ and $G(t,s)$, we have
%
\begin{equation}
\begin{aligned}
% Fu
F^u(t,s,\dot{s}) &=  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot{\color{blue}\underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{u}} :{\color{blue} \underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u})}_{\color{blue}{\boldsymbol{f^u_1}}}} \, d\Omega  + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} \, d\Gamma, \\
% Gu
G^u(t,s) &= 0
\end{aligned}
\end{equation}
%
Note that we have multiple $\vec{f}_0$ functions, each associated with a trial function and an integral over a different domain or boundary.
Each material and boundary condition (except Dirichlet) contribute pointwise functions.
The integral over the domain $\Omega$ is subdivided into integrals over the materials and the integral over the boundary $\Gamma_\tau$ is subdivided into integrals over the Neumann boundaries.
Each bulk constitutive model provides a different pointwise function for the stress tensor $\boldsymbol{\sigma}(\vec{u})$.
With $G=0$ it is clear that we have a formulation that will use implicit time stepping algorithms.

## Jacobian Pointwise Functions

We only have a Jacobian for the LHS:
%
\begin{equation}
\begin{aligned}
J_F^{uu} &= \frac{\partial F^u}{\partial u} = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial u}(-\boldsymbol{\sigma}) \, d\Omega  = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{C} : \frac{1}{2}(\nabla + \nabla^T){\vec{\psi}_\mathit{basis}^{u}}\, d\Omega  = \int_\Omega {\psi_\mathit{trial}^{u}}_{i,k} \, {\color{blue} \underbrace{\color{black}\left( -C_{ikjl} \right)}_{\color{blue}{J_{f3}^{uu}}}} \, {\psi_\mathit{basis}^{u}}_{j,l}\, d\Omega.
\end{aligned}
\end{equation}
