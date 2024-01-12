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
  % LHS displacement
  F^u(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{u}} : {\color{blue}  \underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u},p_f)}_{\color{blue}{\boldsymbol{f}^u_1}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} \, d\Gamma, \\
% RHS displacement
  G^u(t,s) &= 0, \\
% LHS fluid pressure
  F^p(t,s,\dot{s}) &= \int_\Omega  {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\left( \color{black}\frac{\partial \zeta(\vec{u},p_f)}{\partial t} - \gamma(\vec{x},t)\right)}_{\color{blue}{f^p_0}}} + \nabla {\psi_\mathit{trial}^{p}} \cdot {\color{blue}  \underbrace{\color{black}-\vec{q}(p_f)}_{\color{blue}{\vec{f}^p_1}}} \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} ( {\color{blue} \underbrace{\color{black}q_0(\vec{x},t)}_{\color{blue}{f^p_0}}}) \, d\Gamma, \\
% RHS fluid pressure
  G^p(t,s) &= 0, \\
% LHS trace strain
  F^{\epsilon_{v}}(t,s,\dot{s}) &= \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} \cdot {\color{blue}
  \underbrace{\color{black}\left(\nabla \cdot \vec{u} - \epsilon_{v} \right)}_{\color{blue}{f^{\epsilon_{v}}_{0}}}} \, d\Omega. \\
% RHS trace strain
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
  J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_{tshift} \frac{\partial F^u}{\partial \dot{u}} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial u} (- \sigma(\vec{u},p,\epsilon_{v})) \
  d\Omega = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial u} (-(\boldsymbol{C}:\boldsymbol{\varepsilon} -\alpha p \boldsymbol{I})) \ d\Omega \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{C}: \frac{1}{2} (\nabla + \nabla^T) {\vec{\psi}_\mathit{basis}^{u}} \ d\Omega = \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i,k}{\color{blue}  \underbrace{\color{black}\left(-C_{ikjl}\right)}_{\color{blue}{J_{f3}^{uu}}}} {\psi_\mathit{basis}^{u}}_{j,l} \ d\Omega \\
%
% JF_UP
% Jf2up
  J_F^{up} &= \frac{\partial F^u}{\partial p} + s_{tshift} \frac{\partial F^u}{\partial \dot{p}} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial p}(-(\boldsymbol{C}:\boldsymbol{\varepsilon} -\alpha p \boldsymbol{I})) \ d\Omega =
  \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i,j}{\color{blue} \underbrace{\color{black}\left(\alpha \delta_{ij}\right)}_{\color{blue}{J_{f2}^{up}}}} {\psi_\mathit{basis}^{p}} \ d\Omega \\
%
% JF_UE
% Jf2ue
  J_F^{u \epsilon_{v}} &= \frac{\partial F^u}{\partial \epsilon_{v}} + s_{tshift} \frac{\partial F^u}{\partial \dot{\epsilon_{v}}} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial \epsilon_{v}}
  (-\sigma(\vec{u},p,\epsilon_{v})) \ d\Omega = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} :
  \frac{\partial}{\partial \epsilon_{v}} (-(\boldsymbol{C}:\boldsymbol{\varepsilon} -\alpha p \boldsymbol{I})) \ d\Omega \\
  &= \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial \epsilon_{v}} \
  \left[-\left(2 \mu \boldsymbol{\epsilon} + \lambda \boldsymbol{I} \epsilon_{v} - \alpha \boldsymbol{I} p \right) \right] d\Omega =
  \int_{\Omega} {\psi_\mathit{trial}^{u}}_{i,j}{\color{blue} \underbrace{\color{black}\left(-\lambda \delta_{ij} \right)}_{\color{blue}{J_{f2}^{u \epsilon_{v}}}}} {\psi_\mathit{basis}^{\epsilon_{v}}} d\Omega  \\
%
% JF_PU
%
  J_F^{pu} &= \frac{\partial F^p}{\partial u} + s_{tshift} \frac{\partial F^p}{\partial \dot{u}} = 0 \\
%
% JF_PP
% Jf0pp
  J_F^{pp} &= \frac{\partial F^p}{\partial p} + s_{tshift} \frac{\partial F^p}{\partial \dot{p}} =
  \int_{\Omega} \nabla {\psi_\mathit{trial}^{p}} \cdot \frac{\partial}{\partial p} -\left[-\frac{\boldsymbol{k}}{\mu_{f}} \left(\nabla p - \vec{f} \right) \right] \ d\Omega  +
  s_{tshift}\int_{\Omega} {\psi_\mathit{trial}^{p}} \frac{\partial}{\partial \dot{p}} \left[\alpha\dot{\epsilon}_{v} + \frac{\dot{p}}{M} - \gamma\left(\vec{x},t\right)\right] \ d\Omega \\
  &= \int_{\Omega} \nabla \psi_{trial}^ p \left(\frac{\boldsymbol{k}}{\mu_{f}} \nabla \cdot \psi_{basis}^p \right) \ d\Omega +
  \int_{\Omega} {\psi_\mathit{trial}^{p}} \left(s_{tshift} \frac{1}{M}\right) {\psi_\mathit{basis}^{p}} \ d\Omega \\
  &= \int_{\Omega} \psi_{trial,k}^p {\color{blue}\underbrace{\color{black}\left(\frac{\boldsymbol{k}}{\mu_{f}} \delta_{kl}\right)}_{\color{blue}{J_{f3}^{pp}}}} \psi_{basis,l}^p \ d\Omega +\int_{\Omega} {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\color{black}\left(s_{tshift} \frac{1}{M}\right)}_{\color{blue}{J_{f0}^{pp}}}} {\psi_\mathit{basis}^{p}} \ d\Omega \\
%
% JF_PE
% Jf0pe
  J_F^{p\epsilon_{v}} &= \frac{\partial F^p}{\partial \epsilon_{v}} + s_{tshift} \frac{\partial
    F^p}{\partial \dot{\epsilon_{v}}} = \int_{\Omega} {\psi_\mathit{trial}^{p}} {\color{blue}  \underbrace{\color{black}\left(s_{tshift} \alpha \right)}_{\color{blue}{J_{f0}^{p\epsilon_{v}}}}}
    {\psi_\mathit{basis}^{\epsilon_{v}}} \ d\Omega \\
%
% JF_EU
% Jf1eu
  J_F^{\epsilon_{v}u} &= \frac{\partial F^{\epsilon_{v}}}{\partial u} + s_{tshift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{u}} =
  \int_{\Omega} \psi_{trial}^{\epsilon_{v}} \nabla \cdot \vec{\psi}_{basis}^u \ d\Omega = \int_{\Omega}
  {\psi_\mathit{basis}^{\epsilon_{v}}} {\color{blue}  \underbrace{\color{black}\left(\delta_{ij}\right)}_{\color{blue}{J_{f1}^{\epsilon_{v}u}}}}
  {\psi_\mathit{basis}^{u}}_{i,j} \ d\Omega\\
%
% JF_EP
%
  J_F^{\epsilon_{v}p} &= \frac{\partial F^{\epsilon_{v}}}{\partial p} + s_{tshift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{p}} = 0 \\
%
% JF_EE
%
  J_F^{\epsilon_{v}\epsilon_{v}} &= \frac{\partial F^\epsilon_{v}}{\epsilon_{v}} + s_{tshift} \frac{\partial F^{\epsilon_{v}}}{\partial \dot{\epsilon_{v}}} =
  \int_{\Omega} {\psi_\mathit{basis}^{\epsilon_{v}}} {\color{blue}  \underbrace{\color{black}\left(-1\right)}_{\color{blue}{J_{f0}^{\epsilon_{v}\epsilon_{v}}}}} {\psi_\mathit{basis}^{\epsilon_{v}}} \ d\Omega
\end{align}
%

## Porosity State Variable

The default poroelasticity implementation uses a porosity that can vary in space but does not vary in time.
PyLith also includes an implementation with porosity as a state variable that allows it to evolve in time as well as vary in space.

For isothermal conditions {cite:t}`Detournay:Cheng:Poroelasticity:1993` derive a differential equation for the change in porosity:

\begin{equation}
\frac{\partial \phi(\vec{x}, t)}{\partial t} = (\alpha(\vec{x}) - \phi(\vec{x}, t))\left( \dot{\epsilon}_v(\vec{x}, t) + \frac{1-\alpha(\vec{x})}{K_d(\vec{x})} \dot{p}(\vec{x}, t)\right)
\end{equation}

If we use the approximation

\begin{equation}
\frac{\partial \phi(\vec{x}, t)}{\partial t} = \frac{\phi(\vec{x}, t+\Delta t) - \phi(\vec{x}, t)}{\Delta t},
\end{equation}

then we can update the porosity after advancing the solution $\Delta t$ using

\begin{equation}
\phi(\vec{x}, t+\Delta t) = \phi(\vec{x}, t) + \Delta t \left( (\alpha(\vec{x}) - \phi(\vec{x}, t))\left( \dot{\epsilon}_v(\vec{x}, t) + \frac{1-\alpha(\vec{x})}{K_d(\vec{x})} \dot{p}(\vec{x}, t) \right) \right).
\end{equation}

:::{note}
In updating the porosity state variable, we limit values to between 0 and 1, inclusively.
:::

When updating the state variables, PETSc provides the values of the current solution as well as the auxiliary field.
Because the expression for the porosity change depends on the time derivative of the trace strain and fluid pressure, we include those subfields in the solution field.

:::{danger}
In v4.0 we also include the velocity (time derivative of the displacement) in the solution field.
This is not necessary, so we plan to remove it in v5.0.
:::

Our full set of governing equations is

\begin{gather}
% Solution
\vec{s}^T = (\vec{u} \quad p \quad \epsilon_v \quad \vec{v} \quad \dot{p} \quad \dot{\epsilon}_v), \\
% Weak conservation of momentum
  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{\sigma}(\vec{u},p_f) \,
  d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma = 0, \\
% Weak conservation of mass
  \int_\Omega  {\psi_\mathit{trial}^{p}} \left( \frac{\partial \zeta(\vec{u},p_f)}{\partial t} - \gamma(\vec{x},t)\right)
  + \nabla {\psi_\mathit{trial}^{p}} \cdot \left(-\vec{q}(p_f)\right) \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} q_0(\vec{x},t))\, d\Gamma = 0, \text{ and } \\
% Weak vol. strain
  \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} \cdot \left(\nabla \cdot \vec{u} - \epsilon_{v} \right) d\Omega = 0\, \\
% Velocity
\frac{\partial \vec{u}}{\partial t} = \vec{v}\, \\
% Time derivative of fluid pressure
\frac{\partial p}{\partial t} = p_t\, \\
% Time derivative of trace strain
\frac{\partial \epsilon_v}{\partial t} = \epsilon_{vt}.
\end{gather}

### Residual Pointwise Functions

In addition to the residual functions for the case in which porosity does not evolve in time, we have

\begin{align}
% LHS velocity
  F^{v}(t,s,\dot{s}) &= \int_{\Omega} {\vec{\psi}_\mathit{trial}^{v}} \cdot {\color{blue}
  \underbrace{\color{black}\left(\frac{\partial \vec{u}}{\partial t} - \vec{v} \right)}_{\color{blue}{f^{v}_{0}}}} \, d\Omega, \\
% RHS velocity
  G^v(t,s) &= 0, \\
% LHS pressure dot
  F^{p_t}(t,s,\dot{s}) &= \int_{\Omega} {\psi_\mathit{trial}^{p_t}} \cdot {\color{blue}
  \underbrace{\color{black}\left(\frac{\partial p}{\partial t} - p_t \right)}_{\color{blue}{f^{p_t}_{0}}}} \, d\Omega, \\
% RHS pressure dot
  G^{p_t}(t,s) &= 0, \\
% LHS trace strain dot
  F^{\epsilon_{vt}}(t,s,\dot{s}) &= \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{vt}}} \cdot {\color{blue}
  \underbrace{\color{black}\left(\frac{\partial \epsilon_v}{\partial t} - \epsilon_{vt} \right)}_{\color{blue}{f^{\epsilon_{vt}}_{0}}}} \, d\Omega, \\
% RHS trace strain dot
  G^{\epsilon_{vt}}(t,s) &= 0.
\end{align}

### Jacobian Pointwise Functions

The corresponding additional Jacobian functions are

\begin{align}
% JF_vu
  J_F^{vu} &= \frac{\partial F^v}{\partial u} + s_{tshift} \frac{\partial F^v}{\partial \dot{u}} = s_{tshift}, \\
% JF_vv
  J_F^{vv} &= \frac{\partial F^v}{\partial v} + s_{tshift} \frac{\partial F^v}{\partial \dot{v}} = -1, \\
% JF_pdotp
  J_F^{p_t p} &= \frac{\partial F^{p_t}}{\partial p} + s_{tshift} \frac{\partial F^{p_t}}{\partial \dot{p}} = s_{tshift}, \\
% JF_pdotpdot
  J_F^{p_t p_t} &= \frac{\partial F^{p_t}}{\partial {p_t}} + s_{tshift} \frac{\partial F^{p_t}}{\partial \dot{p}_t} = -1, \\
% JF_edpte
  J_F^{\epsilon_{vt} \epsilon_v} &= \frac{\partial F^{\epsilon_{vt}}}{\partial \epsilon_v} + s_{tshift} \frac{\partial F^{\epsilon_{vt}}}{\partial \dot{\epsilon}_t} = s_{tshift}, \\
% JF_edotedot
  J_F^{\epsilon_{vt} \epsilon_{vt}} &= \frac{\partial F^{\epsilon_{vt}}}{\partial {\epsilon_{vt}}} + s_{tshift} \frac{\partial F^{\epsilon_{vt}}}{\partial \dot{\epsilon}_{vt}} = -1, \\
\end{align}
