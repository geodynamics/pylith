# Dynamic

The equation prescribing fault slip is independent of the Lagrange multiplier, so we do not have a system of equations that we can put in
the form $\dot{s} = G^*(t,s)$.
Instead, we have a differential-algebraic set of equations (DAEs), which we solve using an implicit-explicit (IMEX) time integration scheme.
As in the case of dynamic elasticity without faults, we introduce the velocity ($\vec{v}$) as an unknown to turn the elasticity equation into two first order equations.
Our constraint for prescribed slip is
\begin{equation}
  \vec{u}^+ - \vec{u}^- - \vec{d} = \vec{0},
\end{equation}
where $\vec{u}$ is the displacement vector and $\vec{d}$ is the slip vector.
In order to match the order of the time derivative in the elasticity equation, we take the second derivative of the prescribed slip equation with respect to time,
\begin{equation}
  \frac{\partial \vec{v}^+}{\partial t} - \frac{\partial \vec{v}^-}{\partial t} - \frac{\partial^2 \vec{d}}{\partial t^2} = \vec{0}.
\end{equation}
This means that our differential algebraic equations has a differentiation index of 2.

The strong form for our system of equations is:
\begin{gather}
  % Solution
  \vec{s}^T = \left( \vec{u} \quad \vec{v} \quad \vec{\lambda} \right)^T \\
  % Displacement-velocity
  \frac{\partial \vec{u}}{\partial t} = \vec{v}, \\
  % Elasticity
  \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} = \vec{f}(\vec{x},t) + \nabla \cdot \boldsymbol{\sigma}(\vec{u}), \\
  % Neumann BC
  \boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau} \text{ on } \Gamma_\tau. \\
  % Dirichlet BC
  \vec{u} = \vec{u}_0 \text{ on } \Gamma_u, \\
  % Presribed slip
  \frac{\partial \vec{v}^+}{\partial t} - \frac{\partial \vec{v}^-}{\partial t} - \frac{\partial^2 \vec{d}(\vec{x},t)}{\partial t^2} = \vec{0}, \\
  % Momentum balance on fault interface
  \int_{\Omega_f} \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} \, d\Omega = \int_{\Gamma_{f^+}} \boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda} \, d\Gamma + \int_{\Gamma_{f^-}} \boldsymbol{\sigma} \cdot \vec{n} - \vec{\lambda} \, d\Gamma.
\end{gather}
We generate the weak form in the usual way,
\begin{gather}
  % Displacement-velocity
  \int_{\Omega} \vec{\psi}_\mathit{trial}^u \cdot \frac{\partial \vec{u}}{\partial t} \, d\Omega =  \int_{\Omega} \vec{\psi}_\mathit{trial}^u \cdot \vec{v} \, d\Omega, \\
  % Elasticity
  \int_{\Omega} \vec{\psi}_\mathit{trial}^v \cdot \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} \, d\Omega  = \int_\Omega \vec{\psi}_\mathit{trial}^v \cdot \vec{f}(\vec{x},t) + \nabla \vec{\psi}_\mathit{trial}^v : -\boldsymbol{\sigma}(\vec{u}) \, d\Omega + \int_{\Gamma_\tau} \vec{\psi}_\mathit{trial}^v \cdot \vec{\tau}(\vec{x},t) \, d\Gamma \\
  \qquad\qquad + \int_{\Gamma_{f}} \vec{\psi}_\mathit{trial}^{v^+} \cdot \left(-\vec{\lambda}(\vec{x},t)\right) + \vec{\psi}_\mathit{trial}^{v^-} \cdot \left(+\vec{\lambda}(\vec{x},t)\right)\, d\Gamma, \\
  % Prescribed slip
  \int_{\Gamma_f} \vec{\psi}_\mathit{trial}^\lambda \cdot \left(\frac{\partial \vec{v}^+}{\partial t} - \frac{\partial \vec{v}^-}{\partial t} - \frac{\partial^2 \vec{d}(\vec{x},t)}{\partial t^2} \right) \, d\Gamma = 0. \\
  \int_{\Omega_f} \vec{\psi}_\mathit{trial}^\lambda \cdot \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} \, d\Omega = \int_{\Gamma_{f^+}} \vec{\psi}_\mathit{trial}^\lambda \cdot \left( \boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda} \right) \, d\Gamma + \int_{\Gamma_{f^-}} \vec{\psi}_\mathit{trial}^\lambda \cdot \left( \boldsymbol{\sigma} \cdot \vec{n} - \vec{\lambda} \right) \, d\Gamma.
\end{gather}

For compatibility with PETSc TS IMEX implementations, we need $\dot{\vec{s}}$ on the LHS for the explicit part (displacement-velocity and elasticity equations) and we need $\vec{\lambda}$ in the equation for the implicit part (prescribed slip equation).
We first focus on the explicit part and select numerical quadrature that yields a lumped mass matrix, $M$, so that we have

```{math}
:label: eqn:dynamic:prescribed:slip:weak:form
\begin{gathered}
  % Displacement-velocity
  \frac{\partial \vec{u}}{\partial t} = M_u^{-1} \int_{\Omega} \vec{\psi}_\mathit{trial}^u \cdot \vec{v} \, d\Omega, \\
  % Elasticity
  \frac{\partial \vec{v}}{\partial t} = M_v^{-1} \int_\Omega \vec{\psi}_\mathit{trial}^v \cdot \vec{f}(\vec{x},t) + \nabla \vec{\psi}_\mathit{trial}^v : -\boldsymbol{\sigma}(\vec{u}) \, d\Omega + M_v^{-1} \int_{\Gamma_\tau} \vec{\psi}_\mathit{trial}^v \cdot \vec{\tau}(\vec{x},t) \, d\Gamma \\
  \qquad\qquad+ M_{v^+}^{-1} \int_{\Gamma_{f}} \vec{\psi}_\mathit{trial}^{v^+} \cdot \left(-\vec{\lambda}(\vec{x},t)\right) \, d\Gamma + M_{v^-}^{-1} \int_{\Gamma_{f}}\vec{\psi}_\mathit{trial}^{v^-} \cdot \left(+\vec{\lambda}(\vec{x},t)\right) \, d\Gamma, \\
  M_u = \mathit{Lump}\left( \int_\Omega \psi_{\mathit{trial}_i}^u \delta_{ij} \psi_{\mathit{basis}_j}^u \, d\Omega \right), \\
  M_v = \mathit{Lump}\left( \int_\Omega \psi_{\mathit{trial}_i}^v \rho(\vec{x}) \delta_{ij} \psi_{\mathit{basis}_j}^v \, d\Omega \right).
\end{gathered}
```

For the implicit part, we can separate the integration of the weak form for negative and positive sides of the fault interface, which yields
\begin{gather}
  M_{v^+} \frac{\partial \vec{v}^+}{\partial t} = \int_{\Gamma_{f^+}} \vec{\psi}_\mathit{trial}^\lambda \cdot \left( \boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda} \right) \, d\Gamma, \\
  M_{v^-} \frac{\partial \vec{v}^-}{\partial t} = \int_{\Gamma_{f^-}} \vec{\psi}_\mathit{trial}^\lambda \cdot \left( \boldsymbol{\sigma} \cdot \vec{n} - \vec{\lambda} \right) \, d\Gamma.
\end{gather}
Using these equations to substitute in the expressions for the time derivative of the velocity on the negative and positive sides of the fault into the prescribed slip constraint equation yields

```{math}
:label: eqn:elasticity:prescribed:slip:dynamic:DAE:weak:form
  M_{v^+}^{-1} \int_{\Gamma_f^+} \vec{\psi}_\mathit{trial}^\lambda \cdot \left(\boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda}\right) \, d\Gamma + M_{v^-}^{-1} \int_{\Gamma_f^-} \vec{\psi}_\mathit{trial}^\lambda \cdot \left( -\boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda} \right) \, d\Gamma + \int_{\Gamma_f} \vec{\psi}_\mathit{trial}^\lambda \cdot \left(-\frac{\partial^2 \vec{d}}{\partial t^2} \right) \, d\Gamma = \vec{0}.
```

## Residual Pointwise Functions

Combining the explicit parts of the weak form in equation {math:numref}`eqn:dynamic:prescribed:slip:weak:form` with the implicit part of the weak form in equation {math:numref}`eqn:elasticity:prescribed:slip:dynamic:DAE:weak:form` and identifying $F(t,s,\dot{s})$ and $G(t,s)$, we have
\begin{gather}
% Fu
  F^u(t,s,\dot{s}) = \frac{\partial \vec{u}}{\partial t} \\
% Fv
  F^v(t,s,\dot{s}) = \frac{\partial \vec{v}}{\partial t} \\
% Fl
  F^\lambda(t,s,\dot{s}) = {\color{blue}\underbrace{\color{black}M_{v^+}^{-1}}_{\color{blue}{c^+}}} \int_{\Gamma_f^+} \vec{\psi}_\mathit{trial}^\lambda \cdot {\color{blue}\underbrace{\color{black}\left(\boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda}\right)}_{\color{blue}{f^\lambda_0}}} \, d\Gamma + {\color{blue}\underbrace{\color{black}M_{v^-}^{-1}}_{\color{blue}{c^-}}} \int_{\Gamma_f^-} \vec{\psi}_\mathit{trial}^\lambda \cdot {\color{blue}\underbrace{\color{black}\left( -\boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda} \right)}_{\color{blue}{f^\lambda_0}}} \, d\Gamma + \int_{\Gamma_f} \vec{\psi}_\mathit{trial}^\lambda \cdot {\color{blue}\underbrace{\color{black}-\frac{\partial^2 \vec{d}}{\partial t^2}}_{\color{blue}{f^\lambda_0}}} \, d\Gamma \\
% Gu
  G^u(t,s) = {\color{blue}\underbrace{\color{black}M_{u}^{-1}}_{\color{blue}{c}}} \int_\Omega \vec{\psi}_\mathit{trial}^u \cdot {\color{blue}\underbrace{\color{black}\vec{v}}_{\color{blue}{\vec{g}^u_0}}} \, d\Omega, \\
 % Gv
  G^v(t,s) =  {\color{blue}\underbrace{\color{black}M_{v}^{-1}}_{\color{blue}{c}}} \left( \int_\Omega \vec{\psi}_\mathit{trial}^v \cdot {\color{blue}\underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{g}^v_0}}} + \nabla \vec{\psi}_\mathit{trial}^v : {\color{blue}\underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u})}_{\color{blue}{\boldsymbol{g^v_1}}}} \, d\Omega + \int_{\Gamma_\tau} \vec{\psi}_\mathit{trial}^v \cdot {\color{blue}\underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{g}^v_0}}} \, d\Gamma, + \int_{\Gamma_{f}} \vec{\psi}_\mathit{trial}^{v^+} \cdot {\color{blue}\underbrace{\color{black}\left(-\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}{\vec{g}^v_0}}} + \vec{\psi}_\mathit{trial}^{v^-} \cdot {\color{blue}\underbrace{\color{black}\left(+\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}{\vec{g}^v_0}}} \, d\Gamma \right), \\
% Gl
  G^\lambda(t,s) = 0
\end{gather}

The integrals for the explicit part are all weighted by the inverse of the lumped mass matrix.
For the implicit part, only the integrals over the positive and negative sides of the fault are weighted by the inverse of the lumped mass matrix.

## Jacobian Pointwise Functions

For the explicit part we have pointwise functions for computing the lumped LHS Jacobian. These are exactly the same pointwise functions as in the dynamic case without a fault,
\begin{align}
  % J_F uu
  J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{u}} =
             \int_\Omega \psi_{\mathit{trial}^u_i} {\color{blue}\underbrace{\color{black}s_\mathit{tshift} \delta_{ij}}_{\color{blue}{J^{uu}_{f0}}}} \psi_{\mathit{basis}^u_j}  \, d\Omega, \\
  % J_F vv
  J_F^{vv} &= \frac{\partial F^v}{\partial v} + s_\mathit{tshift} \frac{\partial F^v}{\partial \dot{v}} =
             \int_\Omega \psi_{\mathit{trial}^v_i} {\color{blue}\underbrace{\color{black}\rho(\vec{x}) s_\mathit{tshift} \delta_{ij}}_{\color{blue}{J ^{vv}_{f0}}}} \psi_{\mathit{basis}^v_j} \, d\Omega
\end{align}
For the implicit part, we have pointwise functions for the LHS Jacobians associated with the prescribed slip,
\begin{gather}
  % J_F lu
  J_F^{\lambda u} = \frac{\partial F^\lambda}{\partial u} + s_\mathit{tshift} \frac{\partial F^\lambda}{\partial \dot{u}} = {\color{blue}\underbrace{\color{black}M_{v^+}^{-1}}_{\color{blue}{c^+}}} \int_{\Gamma_{f^+}} \psi_{\mathit{trial}_i}^\lambda {\color{blue}\underbrace{\color{black} C_{ikjl} n_k}_{\color{blue}{J^{\lambda u}_{f1}}}} \psi_{\mathit{basis}_{j,l}}^u \, d\Gamma + {\color{blue}\underbrace{\color{black}M_{v^-}^{-1}}_{\color{blue}{c^-}}} \int_{\Gamma_{f^-}} \psi_{\mathit{trial}_i}^\lambda {\color{blue}\underbrace{\color{black}- C_{ikjl} n_k}_{\color{blue}{J^{\lambda u}_{f1}}}} \psi_{\mathit{basis}_{j,l}}^u \, d\Gamma, \\
% J_F ll
  J_F^{\lambda \lambda} = \frac{\partial F^\lambda}{\partial \lambda} + s_\mathit{tshift} \frac{\partial F^\lambda}{\partial \dot{\lambda}} = {\color{blue}\underbrace{\color{black}M_{v^+}^{-1}}_{\color{blue}{c^+}}} \int_{\Gamma_{f^+}} \psi_{\mathit{trial}_i}^\lambda {\color{blue}\underbrace{\color{black} \delta_{ij}}_{\color{blue}{J^{\lambda\lambda}_{f0}}}} \psi_{\mathit{basis}_j}^\lambda \, d\Gamma + {\color{blue}\underbrace{\color{black}M_{v^-}^{-1}}_{\color{blue}{c^-}}} \int_{\Gamma_{f^-}} \psi_{\mathit{trial}_i}^\lambda {\color{blue}\underbrace{\color{black} \delta_{ij}}_{\color{blue}{J^{\lambda\lambda}_{f0}}}} \psi_{\mathit{basis}_j}^\lambda \, d\Gamma
\end{gather}

% End of file
