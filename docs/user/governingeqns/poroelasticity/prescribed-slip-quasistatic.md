# Quasistatic

The strong form of the quasistatic boundary value problem may be expressed as
%
\begin{gather}
% Solution
  \vec{s}^{T} = \left( \vec{u} \quad p \quad \epsilon_v \quad \vec{\lambda} \right), \\
% Elasticity
  \vec{f}(t) + \nabla \cdot \boldsymbol{\sigma}(\vec{u},p) = \vec{0} \text{ in } \Omega, \\
% Pressure
  \frac{\partial \zeta(\vec{u},p)}{\partial t} - \gamma(\vec{x},t) + \nabla \cdot \vec{q}(p) = 0 \text{ in } \Omega, \\
% Vol. Strain
  \nabla \cdot \vec{u} - \epsilon_{v} = 0 \text{ in } \Omega, \\
% Neumann traction
  \boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on } \Gamma_{\tau}, \\
% Fault tractions
  \boldsymbol{\sigma} \cdot \vec{n} = -\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^+}, \\
  \boldsymbol{\sigma} \cdot \vec{n} = +\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^-}, \\
% Neumann flow
  \vec{q} \cdot \vec{n} = q_0(\vec{x}, t) \text{ on } \Gamma_{q}, \\
% Dirichlet displacement
  \vec{u} = \vec{u}_0(\vec{x}, t) \text{ on } \Gamma_{u}, \\
% Dirichlet pressure
  p = p_0(\vec{x},t) \text{ on } \Gamma_{p}, \text{ and} \\
  % Prescribed slip
  -\vec{u}^+ + \vec{u}^- + \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f.
\end{gather}
%
We place all terms for the elasticity, pressure, an volumetric strain equations on the left-hand-side, consistent with PETSc TS implicit time stepping.
We create the weak form by taking the dot product with the trial functions ${\vec{\psi}_\mathit{trial}^{u}}$, ${\psi_\mathit{trial}^{p}}$, ${\psi_\mathit{trial}^{\epsilon_{v}}}$, or $\vec{\psi}_\mathit{trial}^\lambda$ and integrating over the domain:
%
\begin{gather}
% Weak conservation of momentum
  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \left( \vec{f}(\vec{x},t) + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} (\vec{u},p) \right) \, d\Omega = 0, \\
% Weak conservation of mass
  \int_\Omega  {\psi_\mathit{trial}^{p}} \left( \frac{\partial \zeta(\vec{u},p)}{\partial t} - \gamma(\vec{x},t) + \nabla \cdot \vec{q}(p)\right) \, d\Omega = 0,\\
% Weak vol. strain
  \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}}\cdot \left( \nabla \cdot \vec{u} - \epsilon_v \right) \, d\Omega \\
  % Prescribed slip
  \int_{\Gamma_{f}} \vec{\psi}_\mathit{trial}^\lambda \cdot \left(
    -\vec{u}^+ + \vec{u}^- + \vec{d}(\vec{x},t) \right) \, d\Gamma = 0.
\end{gather}
%
Applying the divergence theorem to the first two equations and incorporating the Neumann boundary conditions yields
%
\begin{gather}
% Weak conservation of momentum
  \int_\Omega \vec{\psi_\mathit{trial}^{u}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{\sigma}(\vec{u},p_f) \, d\Omega
  + \int_{\Gamma_\tau} \vec{\psi_\mathit{trial}^{u}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma
  + \int_{\Gamma_{f}} \vec{\psi}_\mathit{trial}^{u^+} \cdot \left(-\vec{\lambda}(\vec{x},t)\right)
  + \vec{\psi}_\mathit{trial}^{u^-} \cdot \left(+\vec{\lambda}(\vec{x},t)\right)\, d\Gamma
   = 0, \\
% Weak conservation of mass
  \int_\Omega  {\psi_\mathit{trial}^{p}} \left( \frac{\partial \zeta(\vec{u},p_f)}{\partial t} - \gamma(\vec{x},t)\right)
  + \nabla {\psi_\mathit{trial}^{p}} \cdot \left(-\vec{q}(p_f)\right) \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} q_0(\vec{x},t))\, d\Gamma = 0, \\
% Weak vol. strain
  \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} \cdot \left(\nabla \cdot \vec{u} - \epsilon_{v} \right) d\Omega = 0, \text{ and } \\
  % Prescribed slip
  \int_{\Gamma_{f}} \vec{\psi}_\mathit{trial}^\lambda \cdot \left(
    -\vec{u}^+ + \vec{u}^- + \vec{d}(\vec{x},t) \right) \, d\Gamma = 0.
\end{gather}
%

## Residual Pointwise Functions

Identifying $F(t,s,\dot{s})$ and $G(t,s)$ we have
%
\begin{align}
  % LHS displacement
  F^u(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{u}} : {\color{blue}  \underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u},p_f)}_{\color{blue}{\boldsymbol{f}^u_1}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} \, d\Gamma + \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{u^+}} \cdot{\color{blue}\underbrace{\color{black}\left(-\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}{\vec{f}^u_0}}} + {\vec{\psi}_\mathit{trial}^{u^-}} \cdot{\color{blue}\underbrace{\color{black}\left(+\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}{\vec{f}^u_0}}}\, d\Gamma, \\
% LHS fluid pressure
  F^p(t,s,\dot{s}) &= \int_\Omega  {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\left( \color{black}\frac{\partial \zeta(\vec{u},p_f)}{\partial t} - \gamma(\vec{x},t)\right)}_{\color{blue}{f^p_0}}} + \nabla {\psi_\mathit{trial}^{p}} \cdot {\color{blue}  \underbrace{\color{black}-\vec{q}(p_f)}_{\color{blue}{\vec{f}^p_1}}} \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} ( {\color{blue} \underbrace{\color{black}q_0(\vec{x},t)}_{\color{blue}{f^p_0}}}) \, d\Gamma, \\
% LHS trace strain
  F^{\epsilon_{v}}(t,s,\dot{s}) &= \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} \cdot {\color{blue}
  \underbrace{\color{black}\left(\nabla \cdot \vec{u} - \epsilon_{v} \right)}_{\color{blue}{f^{\epsilon_{v}}_{0}}}} \, d\Omega. \\
% LHS Lagrange multiplier
F^\lambda(t,s,\dot{s}) &= \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot{\color{blue}\underbrace{\color{black}\left(-\vec{u}^+ + \vec{u}^- + \vec{d}(\vec{x},t) \right)}_{\color{blue}{\vec{f}^\lambda_0}}} \, d\Gamma, \\
% RHS displacement
  G^u(t,s) &= 0, \\
% RHS fluid pressure
  G^p(t,s) &= 0, \\
% RHS trace strain
  G^{\epsilon_v} &= 0, \\
% RHS Lagrange multiplier
  G^{\lambda} &= 0.
\end{align}

## Jacobian Pointwise Functions

Four fields yields potentially 16 Jacobian pointwise functions for the LHS:
%
\begin{align}
%
% Jf3uu
  J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_{tshift} \frac{\partial F^u}{\partial \dot{u}} \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial u} (- \sigma(\vec{u},p,\epsilon_{v})) \, d\Omega \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial u} (-(\boldsymbol{C}:\boldsymbol{\varepsilon} -\alpha p \boldsymbol{I})) \, d\Omega \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{C}: \frac{1}{2} (\nabla + \nabla^T) {\vec{\psi}_\mathit{basis}^{u}} \, d\Omega \\
  &= \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i,k}{\color{blue} \underbrace{\color{orange}\left(-C_{ikjl}\right)}_{\color{blue}{J_{f3}^{uu}}}} {\psi_\mathit{basis}^{u}}_{j,l} \, d\Omega, \\
%
% Jf2up
  J_F^{up} &= \frac{\partial F^u}{\partial p} + s_{tshift} \frac{\partial F^u}{\partial \dot{p}} \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial p}(-(\boldsymbol{C}:\boldsymbol{\varepsilon} -\alpha p \boldsymbol{I})) \, d\Omega \\
  &= \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i,j}{\color{blue} \underbrace{\color{orange}\left(\alpha \delta_{ij}\right)}_{\color{blue}{J_{f2}^{up}}}} {\psi_\mathit{basis}^{p}} \, d\Omega, \\
%
% Jf2ue
  J_F^{u \epsilon_{v}} &= \frac{\partial F^u}{\partial \epsilon_{v}} + s_{tshift} \frac{\partial F^u}{\partial \dot{\epsilon_{v}}} \\
&= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial \epsilon_{v}} (-\sigma(\vec{u},p,\epsilon_{v})) \, d\Omega \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial \epsilon_{v}} (-(\boldsymbol{C}:\boldsymbol{\varepsilon} -\alpha p \boldsymbol{I})) \, d\Omega, \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial \epsilon_{v}} \, \left[-\left(2 \mu \boldsymbol{\epsilon} + \lambda \boldsymbol{I} \epsilon_{v} - \alpha \boldsymbol{I} p \right) \right] d\Omega \\
  &= \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i,j}{\color{blue} \underbrace{\color{orange}\left(-\lambda \delta_{ij} \right)}_{\color{blue}{J_{f2}^{u \epsilon_{v}}}}} {\psi_\mathit{basis}^{\epsilon_{v}}} \, d\Omega,  \\
% J_F ul
J_F^{u\lambda} &= \frac{\partial F^u}{\partial \lambda} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{\lambda}} \\
&= \int_{\Gamma_{f}} {\psi_\mathit{trial}^{u^+}}_i {\color{blue}\underbrace{\color{orange}\left(-\delta_{ij}\right)}_{\color{blue}{J^{u\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j + {\psi_\mathit{trial}^{u^-}}_i {\color{blue}\underbrace{\color{orange}\left(+\delta_{ij}\right)}_{\color{blue}{J^{u\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j \, d\Gamma, \\
%
% Jf_pu
  J_F^{pu} &= \frac{\partial F^p}{\partial u} + s_{tshift} \frac{\partial F^p}{\partial \dot{u}} = 0, \\
%
% Jf0pp
  J_F^{pp} &= \frac{\partial F^p}{\partial p} + s_{tshift} \frac{\partial F^p}{\partial \dot{p}} \\
 &= \int_{\Omega} \nabla {\psi_\mathit{trial}^{p}} \cdot \frac{\partial}{\partial p} -\left[-\frac{\boldsymbol{k}}{\mu_{f}} \left(\nabla p - \vec{f} \right) \right] \, d\Omega + s_{tshift}\int_{\Omega} {\psi_\mathit{trial}^{p}} \frac{\partial}{\partial \dot{p}} \left[\alpha\dot{\epsilon}_{v} + \frac{\dot{p}}{M} - \gamma\left(\vec{x},t\right)\right] \, d\Omega \\
  &= \int_{\Omega} \nabla \psi_{trial}^ p \left(\frac{\boldsymbol{k}}{\mu_{f}} \nabla \cdot \psi_{basis}^p \right) \, d\Omega + \int_{\Omega} {\psi_\mathit{trial}^{p}} \left(s_{tshift} \frac{1}{M}\right) {\psi_\mathit{basis}^{p}} \, d\Omega \\
  &= \int_{\Omega} \psi_{trial,k}^p {\color{blue}\underbrace{\color{orange}\left(\frac{\boldsymbol{k}}{\mu_{f}} \delta_{kl}\right)}_{\color{blue}{J_{f3}^{pp}}}} \psi_{basis,l}^p \ d\Omega +\int_{\Omega} {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\color{orange}\left(s_{tshift} \frac{1}{M}\right)}_{\color{blue}{J_{f0}^{pp}}}} {\psi_\mathit{basis}^{p}} \, d\Omega, \\
%
% Jf0pe
  J_F^{p\epsilon_{v}} &= \frac{\partial F^p}{\partial \epsilon_{v}} + s_{tshift} \frac{\partial F^p}{\partial \dot{\epsilon_{v}}} \\
 &= \int_{\Omega} {\psi_\mathit{trial}^{p}} {\color{blue}  \underbrace{\color{orange}\left(s_{tshift} \alpha \right)}_{\color{blue}{J_{f0}^{p\epsilon_{v}}}}} {\psi_\mathit{basis}^{\epsilon_{v}}} \, d\Omega, \\
%
% Jf_pl
  J_F^{pl} &= \frac{\partial F^p}{\partial \lambda} + s_{tshift} \frac{\partial F^p}{\partial \dot{\lambda}} = 0, \\
%
% Jf1eu
  J_F^{\epsilon_{v}u} &= \frac{\partial F^{\epsilon_{v}}}{\partial u} + s_{tshift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{u}} \\
 &= \int_{\Omega} \psi_{trial}^{\epsilon_{v}} \nabla \cdot \vec{\psi}_{basis}^u \, d\Omega \\
 &= \int_{\Omega} {\psi_\mathit{basis}^{\epsilon_{v}}} {\color{blue}  \underbrace{\color{orange}\left(\delta_{ij}\right)}_{\color{blue}{J_{f1}^{\epsilon_{v}u}}}} {\psi_\mathit{basis}^{u}}_{i,j} \, d\Omega, \\
%
% Jf_ep
  J_F^{\epsilon_{v}p} &= \frac{\partial F^{\epsilon_{v}}}{\partial p} + s_{tshift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{p}} = 0 \\
%
% Jf_ee
  J_F^{\epsilon_{v}\epsilon_{v}} &= \frac{\partial F^\epsilon_{v}}{\epsilon_{v}} + s_{tshift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{\epsilon_{v}}} \\
  &= \int_{\Omega} {\psi_\mathit{basis}^{\epsilon_{v}}} {\color{blue}  \underbrace{\color{orange}\left(-1\right)}_{\color{blue}{J_{f0}^{\epsilon_{v}\epsilon_{v}}}}} {\psi_\mathit{basis}^{\epsilon_{v}}} \, d\Omega, \\
%
% Jf_el
  J_F^{\epsilon_{v}\lambda} &= \frac{\partial F^{\epsilon_{v}}}{\partial \lambda} + s_{tshift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{\lambda}} = 0 \\
%
% J_F lu
J_F^{\lambda u} &= \frac{\partial F^\lambda}{\partial u} + s_\mathit{tshift} \frac{\partial F^\lambda}{\partial \dot{u}} \\
&= \int_{\Gamma_{f}} {\psi_\mathit{trial}^{\lambda}}_i {\color{blue}\underbrace{\color{orange}\left(-\delta_{ij}\right)}_{\color{blue}{J^{\lambda u}_{f0}}}} {\psi_\mathit{basis}^{u^+}}_j  + {\psi_\mathit{trial}^{\lambda}}_i {\color{blue}\underbrace{\color{orange}\left(+\delta_{ij}\right)}_{\color{blue}{J^{\lambda u}_{f0}}}} {\psi_\mathit{basis}^{u^-}}_j \, d\Gamma, \\
%
% Jf_lp
  J_F^{lp} &= \frac{\partial F^l}{\partial p} + s_{tshift} \frac{\partial F^l}{\partial \dot{p}} = 0, \\
%
% Jf_le
  J_F^{l\epsilon_{v}} &= \frac{\partial F^l}{\partial \epsilon_{v}} + s_{tshift} \frac{\partial F^l}{\partial \dot{\epsilon_{v}}} = 0, \\
%
% Jf_ll
  J_F^{ll} &= \frac{\partial F^l}{\partial \lambda} + s_{tshift} \frac{\partial F^l}{\partial \dot{\lambda}} = 0, \\
\end{align}
