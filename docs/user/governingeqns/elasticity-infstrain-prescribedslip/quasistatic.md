# Quasistatic

As in the case of elasticity without faults, we first consider the quasistatic case in which we neglect the inertial term ($\rho \frac{\partial \vec{v}}{\partial t} \approx \vec{0}$).
We place all of the terms in the elasticity equation on the LHS, consistent with implicit time stepping.
Our equation of the conservation of momentum on the fault interface reduces to
\begin{equation}
  \int_{\Gamma_{f^+}} \boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda} \, d\Gamma + \int_{\Gamma_{f^-}} \boldsymbol{\sigma} \cdot \vec{n} - \vec{\lambda} \, d\Gamma = 0.
\end{equation}
We enforce this equation on each portion of the fault interface along with our prescribed slip constraint, which leads to
\begin{gather}
  \boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda} = \vec{0} \text{ on } \Gamma_{f^+}, \\
  \boldsymbol{\sigma} \cdot \vec{n} - \vec{\lambda} = \vec{0}\text{ on } \Gamma_{f^-}, \\
  \vec{u}^+ - \vec{u}^- - \vec{d}(\vec{x},t) = \vec{0},  
\end{gather}

Our solution vector consists of both displacements and Lagrange multipliers, and the strong form for the system of equations is
\begin{gather}
  % Solution
  \vec{s}^T = \left( \vec{u} \quad \vec{\lambda} \right)^T \\
  % Elasticity
  \vec{f}(\vec{x},t) + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u}) = \vec{0} \text{ in }\Omega, \\
  % Neumann
  \boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on }\Gamma_\tau, \\
  % Dirichlet
  \vec{u} = \vec{u}_0(\vec{x},t) \text{ on }\Gamma_u, \\
  % Prescribed slip
  \vec{u}^+ - \vec{u}^- - \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f,  \\
  \boldsymbol{\sigma} \cdot \vec{n} = -\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^+}, \text{ and}\\
  \boldsymbol{\sigma} \cdot \vec{n} = +\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^-}.
\end{gather}
We create the weak form by taking the dot product with the trial function $\vec{\psi}_\mathit{trial}^u$ or $\vec{\psi}_\mathit{trial}^\lambda$ and integrating over the domain.
After using the divergence theorem and incorporating the Neumann boundary and fault interface conditions, we have
\begin{gather}
  % Elasticity
  \int_\Omega \vec{\psi}_\mathit{trial}^u \cdot \vec{f}(\vec{x},t) + \nabla \vec{\psi}_\mathit{trial}^u : -\boldsymbol{\sigma}(\vec{u}) \, d\Omega
  + \int_{\Gamma_\tau} \vec{\psi}_\mathit{trial}^u \cdot \vec{\tau}(\vec{x},t) \, d\Gamma,
  + \int_{\Gamma_{f}} \vec{\psi}_\mathit{trial}^{u^+} \cdot \left(-\vec{\lambda}(\vec{x},t)\right)
  + \vec{\psi}_\mathit{trial}^{u^-} \cdot \left(+\vec{\lambda}(\vec{x},t)\right)\, d\Gamma = 0\\
  % Prescribed slip
  \int_{\Gamma_{f}} \vec{\psi}_\mathit{trial}^\lambda \cdot \left(
    -\vec{u}^+ + \vec{u}^- + \vec{d}(\vec{x},t) \right) \, d\Gamma = 0.
\end{gather}
We have chosen the signs in the last equation so that the Jacobian will be symmetric with respect to the Lagrange multiplier.
We solve the system of equations using implicit time stepping, making use of residuals functions and Jacobians for the LHS.

## Residual Pointwise Functions

Identifying $F(t,s,\dot{s})$ and $G(t,s)$, we have
%
\begin{equation}
\begin{aligned}
% Fu
F^u(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot{\color{blue}\underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{u}} : {\color{blue}\underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u})}_{\color{blue}{\boldsymbol{f^u_1}}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot{\color{blue}\underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} \, d\Gamma + \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{u^+}} \cdot{\color{blue}\underbrace{\color{black}\left(-\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}{\vec{f}^u_0}}} + {\vec{\psi}_\mathit{trial}^{u^-}} \cdot{\color{blue}\underbrace{\color{black}\left(+\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}{\vec{f}^u_0}}}\, d\Gamma \\
% Fl
F^\lambda(t,s,\dot{s}) &= \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot{\color{blue}\underbrace{\color{black}\left(-\vec{u}^+ + \vec{u}^- + \vec{d}(\vec{x},t) \right)}_{\color{blue}{\vec{f}^\lambda_0}}} \, d\Gamma, \\
% Gu
G^u(t,s) &= 0 \\
% Gl
G^\lambda(t,s) &= 0
\end{aligned}
\end{equation}
%
Compared to the quasistatic elasticity case without a fault, we have simply added additional pointwise functions associated with the fault.
Our fault implementation does not change the formulation for the materials or external Dirichlet or Neumann boundary conditions.

## Jacobian Pointwise Functions

The LHS Jacobians are:
%
\begin{equation}
\begin{aligned}
% J_F uu
J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{u}} = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{C} : \frac{1}{2}(\nabla + \nabla^T){\vec{\psi}_\mathit{basis}^{u}} \, d\Omega = \int_\Omega {\psi_\mathit{trial}^{u}}_{i,k} \,{\color{blue}\underbrace{\color{black}\left( -C_{ikjl} \right)}_{\color{blue}{J_{f3}^{uu}}}} \, {\psi_\mathit{basis}^{u}}_{j,l}\, d\Omega \\
% J_F ul
J_F^{u\lambda} &= \frac{\partial F^u}{\partial \lambda} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{\lambda}} = \int_{\Gamma_{f}} {\psi_\mathit{trial}^{u^+}}_i {\color{blue}\underbrace{\color{black}\left(-\delta_{ij}\right)}_{\color{blue}{J^{u\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j + {\psi_\mathit{trial}^{u^-}}_i {\color{blue}\underbrace{\color{black}\left(+\delta_{ij}\right)}_{\color{blue}{J^{u\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j\, d\Gamma, \\
% J_F lu
J_F^{\lambda u} &= \frac{\partial F^\lambda}{\partial u} + s_\mathit{tshift} \frac{\partial F^\lambda}{\partial \dot{u}} = \int_{\Gamma_{f}} {\psi_\mathit{trial}^{\lambda}}_i {\color{blue}\underbrace{\color{black}\left(-\delta_{ij}\right)}_{\color{blue}{J^{\lambda u}_{f0}}}} {\psi_\mathit{basis}^{u^+}}_j  + {\psi_\mathit{trial}^{\lambda}}_i {\color{blue}\underbrace{\color{black}\left(+\delta_{ij}\right)}_{\color{blue}{J^{\lambda u}_{f0}}}} {\psi_\mathit{basis}^{u^-}}_j \, d\Gamma, \\
% J_F ll
J_F^{\lambda \lambda} &= \boldsymbol{0}
\end{aligned}
\end{equation}
%
This LHS Jacobian has the structure
%
\begin{equation}
J_F = \left( \begin{array} {cc} J_F^{uu} & J_F^{u\lambda} \\ J_F^{\lambda u} & 0 \end{array} \right) = \left( \begin{array} {cc} J_F^{uu} & C^T \\ C & 0 \end{array} \right),
\end{equation}  
%
where $C$ contains entries of $\pm 1$ for degrees of freedom on the two sides of the fault.
The Schur complement of $J$ with respect to $J_F^{uu}$ is $-C\left(J_F^{uu}\right)^{-1}C^T$.
