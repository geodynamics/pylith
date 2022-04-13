# Quasistatic

For ease of solution in the quasistatic case, we introduce a third variable in the form of volumetric strain ($\epsilon_v$).
The strong form of the problem may be expressed as
%
\begin{gather}
% Solution
  \vec{s}^{T} = \left(\vec{u} \quad p \quad \epsilon_v\right), \\
% Elasticity
  \vec{f}(t) + \nabla \cdot \boldsymbol{\sigma}(\vec{u},p) = \vec{0} \text{ in } \Omega, \\
% Pressure
  \frac{\partial \zeta(\vec{u},p)}{\partial t} - \gamma(\vec{x},t) + \nabla \cdot \vec{q}(p) = 0 \text{ in } \Omega, \\
% Vol. Strain
  \nabla \cdot \vec{u} - \epsilon_{v} = 0 \text{ in } \Omega, \\
% Neumann traction
  \boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on } \Gamma_{\tau}, \\
% Dirichlet displacement
  \vec{u} = \vec{u}_0(\vec{x}, t) \text{ on } \Gamma_{u}, \\
% Neumann flow
  \vec{q} \cdot \vec{n} = q_0(\vec{x}, t) \text{ on } \Gamma_{q}, \text{ and } \\
% Dirichlet pressure
  p = p_0(\vec{x},t) \text{ on } \Gamma_{p}.
\end{gather}
%
We place all terms for the elasticity, pressure, an volumetric strain equations on the left-hand-side, consistent with PETSc TS implicit time stepping.

We create the weak form by taking the dot product with the trial functions ${\vec{\psi}_\mathit{trial}^{u}}$, ${\psi_\mathit{trial}^{p}}$, and ${\psi_\mathit{trial}^{\epsilon_{v}}}$ and integrating over the domain:
%
\begin{gather}
% Weak conservation of momentum
  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \left( \vec{f}(\vec{x},t) + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} (\vec{u},p) \right) \, d\Omega = 0, \\
% Weak conservation of mass
  \int_\Omega  {\psi_\mathit{trial}^{p}} \left( \frac{\partial \zeta(\vec{u},p)}{\partial t} - \gamma(\vec{x},t) + \nabla \cdot \vec{q}(p)\right) \, d\Omega = 0,\\
% Weak vol. strain
  \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}}\cdot \left( \nabla \cdot \vec{u} - \epsilon_v \right) \, d\Omega.
\end{gather}
%
Applying the divergence theorem to the first two equations and incorporating the Neumann boundary conditions yields
%
\begin{gather}
% Weak conservation of momentum
  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{\sigma}(\vec{u},p_f) \,
  d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma = 0, \\
% Weak conservation of mass
  \int_\Omega  {\psi_\mathit{trial}^{p}} \left( \frac{\partial \zeta(\vec{u},p_f)}{\partial t} - \gamma(\vec{x},t)\right)
  + \nabla {\psi_\mathit{trial}^{p}} \cdot \left(-\vec{q}(p_f)\right) \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} q_0(\vec{x},t))\, d\Gamma = 0, \text{ and } \\
% Weak vol. strain
  \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} \cdot \left(\nabla \cdot \vec{u} - \epsilon_{v} \right) d\Omega = 0
\end{gather}
%

## Residual Pointwise Functions

Identifying $F(t,s,\dot{s})$ and $G(t,s)$ we have
%
\begin{align}
  % Displacement
  F^u(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{u}} :  {\color{blue}  \underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u},p_f)}_{\color{blue}{\boldsymbol{f}^u_1}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} \, d\Gamma, \\
% Pressure
  F^p(t,s,\dot{s}) &= \int_\Omega  {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\left( \color{black}\frac{\partial \zeta(\vec{u},p_f)}{\partial t} - \gamma(\vec{x},t)\right)}_{\color{blue}{f^p_0}}} + \nabla {\psi_\mathit{trial}^{p}} \cdot {\color{blue}  \underbrace{\color{black}-\vec{q}(p_f)}_{\color{blue}{\vec{f}^p_1}}} \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} ( {\color{blue} \underbrace{\color{black}q_0(\vec{x},t)}_{\color{blue}{f^p_0}}}) \, d\Gamma, \\
% Volumetric Strain
  F^{\epsilon_{v}}(t,s,\dot{s}) &= \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} \cdot {\color{blue}
  \underbrace{\color{black}\left(\nabla \cdot \vec{u} - \epsilon_{v} \right)}_{\color{blue}{f^{\epsilon_{v}}_{0}}}} \, d\Omega. \\
  G^u(t,s) &= 0, \\
  G^p(t,s) &= 0, \\
  G^{\epsilon_v} &= 0.
\end{align}
%
## Jacobian Pointwise Functions

Three field results in a potential nine Jacobian pointwise functions for the LHS:
%
\begin{align}
%
% JF_UU
% Jf3uu
  J_F^{uu} &= \frac{\partial F^u}{\partial u} + t_{shift} \frac{\partial F^u}{\partial \dot{u}} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial u} (- \sigma(\vec{u},p,\epsilon_{v})) \
  d\Omega = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial u} (-(\boldsymbol{C}:\boldsymbol{\varepsilon} -\alpha p \boldsymbol{I})) \ d\Omega \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{C}: \frac{1}{2} (\nabla + \nabla^T) {\vec{\psi}_\mathit{basis}^{u}} \ d\Omega = \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i,k}{\color{blue}  \underbrace{\color{black}\left(-C_{ikjl}\right)}_{\color{blue}{J_{f3}^{uu}}}} {\psi_\mathit{basis}^{u}}_{j,l} \ d\Omega \\
%
% JF_UP
% Jf2up
  J_F^{up} &= \frac{\partial F^u}{\partial p} + t_{shift} \frac{\partial F^u}{\partial \dot{p}} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial p}(-(\boldsymbol{C}:\boldsymbol{\varepsilon} -\alpha p \boldsymbol{I})) \ d\Omega =
  \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i,j}{\color{blue} \underbrace{\color{black}\left(\alpha \delta_{ij}\right)}_{\color{blue}{J_{f2}^{up}}}} {\psi_\mathit{basis}^{p}} \ d\Omega \\
%
% JF_UE
% Jf2ue
  J_F^{u \epsilon_{v}} &= \frac{\partial F^u}{\partial \epsilon_{v}} + t_{shift} \frac{\partial F^u}{\partial \dot{\epsilon_{v}}} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial \epsilon_{v}}
  (-\sigma(\vec{u},p,\epsilon_{v})) \ d\Omega = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} :
  \frac{\partial}{\partial \epsilon_{v}} (-(\boldsymbol{C}:\boldsymbol{\varepsilon} -\alpha p \boldsymbol{I})) \ d\Omega \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial \epsilon_{v}} \
  \left[-\left(2 \mu \boldsymbol{\epsilon} + \lambda \boldsymbol{I} \epsilon_{v} - \alpha \boldsymbol{I} p \right) \right] d\Omega =
  \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i,j}{\color{blue} \underbrace{\color{black}\left(-\lambda \delta_{ij} \right)}_{\color{blue}{J_{f2}^{u \epsilon_{v}}}}} {\psi_\mathit{basis}^{\epsilon_{v}}} d\Omega  \\
%
% JF_PU
%
  J_F^{pu} &= \frac{\partial F^p}{\partial u} + t_{shift} \frac{\partial F^p}{\partial \dot{u}} = 0 \\
%
% JF_PP
% Jf0pp
  J_F^{pp} &= \frac{\partial F^p}{\partial p} + t_{shift} \frac{\partial F^p}{\partial \dot{p}} =
  \int_{\Omega} \nabla {\psi_\mathit{trial}^{p}} \cdot \frac{\partial}{\partial p} -\left[-\frac{\boldsymbol{k}}{\mu_{f}} \left(\nabla p - \vec{f} \right) \right] \ d\Omega  +
  t_{shift}\int_{\Omega} {\psi_\mathit{trial}^{p}} \frac{\partial}{\partial \dot{p}} \left[\alpha\dot{\epsilon}_{v} + \frac{\dot{p}}{M} - \gamma\left(\vec{x},t\right)\right] \ d\Omega \\
  &= \int_{\Omega} \nabla \psi_{trial}^ p \left(\frac{\boldsymbol{k}}{\mu_{f}} \nabla \cdot \psi_{basis}^p \right) \ d\Omega +
  \int_{\Omega} {\psi_\mathit{trial}^{p}} \left(t_{shift} \frac{1}{M}\right) {\psi_\mathit{basis}^{p}} \ d\Omega \\
  &= \int_{\Omega} \psi_{trial,k}^p {\color{blue}\underbrace{\color{black}\left(\frac{\boldsymbol{k}}{\mu_{f}} \delta_{kl}\right)}_{\color{blue}{J_{f3}^{pp}}}} \psi_{basis,l}^p \ d\Omega +\int_{\Omega} {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\color{black}\left(t_{shift} \frac{1}{M}\right)}_{\color{blue}{J_{f0}^{pp}}}} {\psi_\mathit{basis}^{p}} \ d\Omega \\
%
% JF_PE
% Jf0pe
  J_F^{p\epsilon_{v}} &= \frac{\partial F^p}{\partial \epsilon_{v}} + t_{shift} \frac{\partial
    F^p}{\partial \dot{\epsilon_{v}}} = \int_{\Omega} {\psi_\mathit{trial}^{p}} {\color{blue}  \underbrace{\color{black}\left(t_{shift} \alpha \right)}_{\color{blue}{J_{f0}^{p\epsilon_{v}}}}}
    {\psi_\mathit{basis}^{\epsilon_{v}}} \ d\Omega \\
%
% JF_EU
% Jf1eu
  J_F^{\epsilon_{v}u} &= \frac{\partial F^{\epsilon_{v}}}{\partial u} + t_{shift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{u}} =
  \int_{\Omega} \psi_{trial}^{\epsilon_{v}} \nabla \cdot \vec{\psi}_{basis}^u \ d\Omega = \int_{\Omega}
  {\psi_\mathit{basis}^{\epsilon_{v}}} {\color{blue}  \underbrace{\color{black}\left(\delta_{ij}\right)}_{\color{blue}{J_{f1}^{\epsilon_{v}u}}}}
  {\psi_\mathit{basis}^{u}}_{i,j} \ d\Omega\\
%
% JF_EP
%
  J_F^{\epsilon_{v}p} &= \frac{\partial F^{\epsilon_{v}}}{\partial p} + t_{shift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{p}} = 0 \\
%
% JF_EE
%
  J_F^{\epsilon_{v}\epsilon_{v}} &= \frac{\partial F^\epsilon_{v}}{\epsilon_{v}} + t_{shift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{\epsilon_{v}}} =
  \int_{\Omega} {\psi_\mathit{basis}^{\epsilon_{v}}} {\color{blue}  \underbrace{\color{black}\left(-1\right)}_{\color{blue}{J_{f0}^{\epsilon_{v}\epsilon_{v}}}}} {\psi_\mathit{basis}^{\epsilon_{v}}} \ d\Omega
\end{align}
%
