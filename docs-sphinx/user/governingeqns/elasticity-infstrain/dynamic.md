# Dynamic

For compatibility with PETSc TS algorithms, we want to turn the second order equation {math:numref}`eqn:elasticity:strong:form` into two first order
equations.
We introduce the velocity as a unknown, $\vec{v}=\frac{\partial u}{\partial t}$, which leads to
%
\begin{equation}
\begin{aligned}
% Displacement-velocity
\frac{\partial \vec{u}}{\partial t} &= \vec{v} \text{ in } \Omega, \\
% Elasticity
\rho(\vec{x}) \frac{\partial\vec{v}}{\partial t} &= \vec{f}(\vec{x},t) + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u}) \text{ in } \Omega, \\
% Neumann
\boldsymbol{\sigma} \cdot \vec{n} &= \vec{\tau}(\vec{x},t) \text{ on } \Gamma_\tau, \\
% Dirichlet
\vec{u} &= \vec{u}_0(\vec{x},t) \text{ on } \Gamma_u.
\end{aligned}
\end{equation}
%
We create the weak form by taking the dot product with the trial function ${\vec{\psi}_\mathit{trial}^{u}}$ or ${\vec{\psi}_\mathit{trial}^{v}}$ and integrating over the domain:
%
\begin{gather}
% Displacement-velocity
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \frac{\partial \vec{u}}{\partial t} \, d\Omega = \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{v} \, d\Omega, \\
% Elasticity
\int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} \, d\Omega = \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \left( \vec{f}(t) + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} (\vec{u}) \right) \, d\Omega.
\end{gather}
%
Using the divergence theorem and incorporating the Neumann boundaries, we can rewrite the second equation as
%
\begin{equation}
\int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} \, d\Omega= \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\boldsymbol{\sigma}(\vec{u}) \, d\Omega  + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma.
\end{equation}
%
For explicit time stepping, we want $F(t,s,\dot{s})=\dot{s}$, so we solve an augmented system in which we multiply the RHS residual function by the inversion of the lumped LHS Jacobian,
%
\begin{gather}
F^*(t,s,\dot{s}) = G^*(t,s) \text{, where} \\
F^*(t,s,\dot{s}) = \dot{s} \text{ and} \\
G^*(t,s) = J_F^{-1} G(t,s).
\end{gather}
%
With the augmented system, we have
%
\begin{gather}
% Displacement-velocity
\frac{\partial \vec{u}}{\partial t}  = M_u^{-1} \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{v} \, d\Omega, \\
% Elasticity
\frac{\partial \vec{v}}{\partial t} = M_v^{-1} \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \left( \vec{f}(t) + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} (\vec{u}) \right) \, d\Omega, \\
% Mu
M_u = \mathit{Lump}\left( \int_\Omega {\psi_\mathit{trial}^{u}}_i \delta_{ij} {\psi_\mathit{basis}^{u}}_j \, d\Omega \right), \\
% Mv
M_v = \mathit{Lump}\left( \int_\Omega {\psi_\mathit{trial}^{v}}_i \rho(\vec{x}) \delta_{ij} {\psi_\mathit{basis}^{v}}_j \, d\Omega \right).
\end{gather}
%
## Residual Pointwise Functions

With explicit time stepping the PETSc TS assumes the LHS is $\dot{s}$, so we only need the RHS residual functions:
%
\begin{equation}
\begin{aligned}
% Gu
G^u(t,s) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}\underbrace{\color{black}\vec{v}}_{\color{blue}{\vec{g}^u_0}}} \, d\Omega, \\
% Gv
G^v(t,s) &=  \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot {\color{blue}\underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{g}^v_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{v}} : {\color{blue}\underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u})}_{\color{blue}{\boldsymbol{g^v_1}}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot{\color{blue}\underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{g}^v_0}}} \, d\Gamma,
\end{aligned}
\end{equation}
%
In the second equation these are the same pointwise functions as in the quasistatic case, only now they are on the RHS instead of the LHS.

## Jacobian Pointwise Functions

These are the pointwise functions associated with $M_u$ and $M_v$ for computing the lumped LHS Jacobian.
We premultiply the RHS residual function by the inverse of the lumped LHS Jacobian while $s_\mathit{tshift}$ remains on the LHS with $\dot{s}$. As a result, we use LHS Jacobian pointwise functions, but set $s_\mathit{tshift}=1$.
The LHS Jacobians are:
%
\begin{equation}
\begin{aligned}
% J_F uu
M_u = J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{u}} = \int_\Omega {\psi_\mathit{trial}^{u}}_i{\color{blue}\underbrace{\color{black}s_\mathit{tshift} \delta_{ij}}_{\color{blue}{J^{uu}_{f0}}}} {\psi_\mathit{basis}^{u}}_j  \, d\Omega, \\
% J_F vv
M_v = J_F^{vv} &= \frac{\partial F^v}{\partial v} + s_\mathit{tshift} \frac{\partial F^v}{\partial \dot{v}} = \int_\Omega {\psi_\mathit{trial}^{v}}_i {\color{blue}  \underbrace{\color{black}\rho(\vec{x}) s_\mathit{tshift} \delta_{ij}}_{\color{blue}{J ^{vv}_{f0}}}} {\psi_\mathit{basis}^{v}}_j \, d\Omega
\end{aligned}
\end{equation}
