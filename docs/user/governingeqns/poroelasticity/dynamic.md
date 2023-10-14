# Dynamic

:::{admonition} TODO @rwalkerlewis
This section will need to be updated for consistency with the dynamic prescribed slip formulation.
:::

For compatibility with PETSc TS algorithms, we want to turn the second order elasticity equation into two first order equations.
We introduce velocity as a unknown, $\vec{v}=\frac{\partial u}{\partial t}$, which leads to a slightly different three field problem,
%
\begin{gather}
% Solution
\vec{s}^{T} = \left(\vec{u} \quad p \quad \vec{v}\right) \\
% Displacement
\frac{\partial \vec{u}}{\partial t} = \vec{v} \text{ in } \Omega \\
% Pressure
\frac{\partial \zeta(\vec{u},p)}{\partial t } - \gamma(\vec{x},t) + \nabla \cdot \vec{q}(p) = 0 \text{ in } \Omega \\
% Velocity
\rho_{b} \frac{\partial \vec{v}}{\partial t} = \vec{f}(\vec{x},t) + \nabla \cdot \boldsymbol{\sigma}(\vec{u},p) \text{ in } \Omega \\
% Neumann traction
\boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on } \Gamma_{\tau} \\
% Dirichlet displacement
\vec{u} = \vec{u}_{0}(\vec{x}, t) \text{ on } \Gamma_{u} \\
% Neumann flow
\vec{q} \cdot \vec{n} = q_{0}(\vec{x}, t) \text{ on } \Gamma_{q} \\
% Dirichlet pressure
p = p_{0}(\vec{x},t) \text{ on } \Gamma_{p}
\end{gather}
%
For compatibility with PETSc TS explicit time stepping algorithms, we need the left hand side to be $F = (t,s,\dot{s}) = \dot{s}$.
We replace the variation of fluid content variable, $\zeta$, with its definition in the conservation of fluid mass equation and solve for the rate of change of pressure,
%
\begin{gather}
  \frac{\partial}{\partial t}\left(\alpha \epsilon_{v} + \frac{p}{M}\right) - \gamma\left(\vec{x},t\right) + \nabla \cdot \vec{q} = 0 \\
  \alpha \dot{\epsilon}_{v} + \frac{\dot{p}}{M} - \gamma \left(\vec{x},t\right) + \nabla \cdot \vec{q} = 0 \\
  \frac{\dot{p}}{M} = \gamma \left(\vec{x},t \right) - \alpha \dot{\epsilon}_{v} -\nabla \cdot \vec{q} \\
  \frac{\dot{p}}{M} = \gamma \left(\vec{x},t \right) - \alpha \left( \nabla \cdot \dot{\vec{u}} \right) -\nabla \cdot \vec{q}.
\end{gather}
%
We write the volumetric strain in terms of displacement, because this dynamic formulation does not include the volumetric strain as an unknown. Note that for poroelastodynamics we use the generalized Darcy's law {cite}`ding2013fundamental` as
%
\begin{equation}
  \vec{q}(p) = -\frac{\boldsymbol{k}}{\mu_{f}}(\nabla p + \rho_{f} \frac{\partial \vec{v}}{\partial t} - \vec{f}_f).
\end{equation}
%
Using trial functions ${\vec{\psi}_\mathit{trial}^{u}}$, ${\psi_\mathit{trial}^{p}}$, and ${\vec{\psi}_\mathit{trial}^{v}}$, and incorporating the Neumann boundary conditions, the weak form may be written as:
%
\begin{align}
  % Displacement
  \int_{\Omega} {\vec{\psi}_\mathit{trial}^{u}} \cdot \left( \frac{\partial \vec{u}}{\partial t} \right)d \Omega &= \int_{\Omega} {\vec{\psi}_\mathit{trial}^{u}} \cdot \left( \vec{v} \right) d \Omega \\
  % Pressure
  \int_{\Omega} {\psi_\mathit{trial}^{p}} \left( \frac{1}{M}\frac{\partial p}{\partial t} \right) d\Omega &=
  \int_{\Omega} {\psi_\mathit{trial}^{p}} \left[\gamma(\vec{x},t) - \alpha \left(\nabla \cdot \dot{\vec{u}}\right) \right]  + \nabla {\psi_\mathit{trial}^{p}} \cdot \vec{q}(p) \ d\Omega +
  \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} (-q_0(\vec{x},t)) \ d\Gamma, \\
  % Velocity
  \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \left( \rho_{b} \frac{\partial
  \vec{v}}{\partial t} \right) \,
  d\Omega &= \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{v}} :
  -\boldsymbol{\sigma} (\vec{u},p_f) \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}}
  \cdot \vec{\tau}(\vec{x},t) \, d\Gamma.
\end{align}
%
## Residual Pointwise Functions

With explicit time stepping the PETSc TS assumes the LHS is $\dot{s}$ , so we only need the RHS residual functions:
%
\begin{align}
% Displacement
  G^u(t,s) &= \int_{\Omega} {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue} \underbrace{\color{black}\vec{v}}_{\color{blue}{\vec{g}_0^u}}} d \Omega, \\
% Pressure
  G^p(t,s) &= \int_\Omega {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\color{black}\left(\gamma(\vec{x},t) - \alpha (\nabla \cdot \dot{\vec{u}})\right)}_{\color{blue}{g^p_0}}} + \nabla {\psi_\mathit{trial}^{p}} \cdot {\color{blue} \underbrace{\color{black}\vec{q}(p_f)}_{\color{blue}{\vec{g}^p_1}}} \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} ({\color{blue} \underbrace{\color{black}-q_0(\vec{x},t)}_{\color{blue}{g^p_0}}}) \, d\Gamma, \\
% Velocity
 G^v(t,s) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot {\color{blue}  \underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{g}^v_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{v}} : {\color{blue}  \underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u},p_f)}_{\color{blue}{\boldsymbol{g}^v_1}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot   {\color{blue}  \underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{g}^v_0}}} \, d\Gamma.
\end{align}
%

## Jacobians Pointwise Functions

These are the pointwise functions associated with $M_{u}$, $M_{p}$, and $M_{v}$ for computing the lumped LHS Jacobian.
We premultiply the RHS residual function by the inverse of the lumped LHS Jacobian while $s_\mathit{tshift}$ remains on the LHS with $\dot{s}$. As a result, we use LHS Jacobian pointwise functions, but set $s_\mathit{tshift} = 1$.
The LHS Jacobians are:
%
\begin{align}
% Displacement
  M_{u} &= J_F^{uu} = \frac{\partial F^u}{\partial u} + s_{tshift} \frac{\partial F^u}{\partial \dot{u}} =
  \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i} {\color{blue}  \underbrace{\color{black}s_{tshift} \delta_{ij}}_{\color{blue}{J^{uu}_{f0}}}} {\psi_\mathit{basis}^{u}}_{j} \, d \Omega \\
% Pressure
  M_{p} &= J_F^{pp} = \frac{\partial F^p}{\partial p} + s_{tshift} \frac{\partial F^p}{\partial \dot{p}} =
  \int_{\Omega} {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\color{black}\left(s_{tshift} \frac{1}{M}\right)}_{\color{blue}{J_{f0}^{pp}}}} {\psi_\mathit{basis}^{p}} \ d\Omega \\
% Velocity
  M_{v} &= J_F^{vv} = \frac{\partial F^v}{\partial v} + s_{tshift} \frac{\partial F^v}{\partial \dot{v}} =
  \int_{\Omega} {\psi_\mathit{trial}^{v}}_{i}{\color{blue}  \underbrace{\color{black}\rho_{b}(\vec{x}) s_{tshift} \delta_{ij}}_{\color{blue}{J^{vv}_{f0}}}} {\psi_\mathit{basis}^{v}}_{j} \  d \Omega
\end{align}
