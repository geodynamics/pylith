# Dynamic

The equation prescribing fault slip is independent of the Lagrange multiplier, so we do not have a system of equations that we can put in the form $\dot{s} = G^*(t,s)$.
Instead, we have a differential-algebraic set of equations (DAEs), which we solve using an implicit-explicit (IMEX) time integration scheme.
As in the case of dynamic elasticity without faults, we introduce the velocity ($\vec{v}$) as an unknown to turn the elasticity equation into two first order equations.
The strong form for our system of equations is:
%
\begin{gather}
% Solution
\vec{s}^T = \left( \vec{u} \quad \vec{v} \quad \vec{\lambda} \right)^T \\
% Displacement-velocity
\frac{\partial \vec{u}}{\partial t} = \vec{v}, \\
% Elasticity
\rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} =
\vec{f}(\vec{x},t) + \nabla \cdot \boldsymbol{\sigma}(\vec{u}), \\
% Neumann BC
\boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau} \text{ on } \Gamma_\tau. \\
% Dirichlet BC
\vec{u} = \vec{u}_0 \text{ on } \Gamma_u, \\
% Presribed slip
\vec{u}^+(\vec{x},t) - \vec{u}^-(\vec{x},t) - \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f,  \\
\boldsymbol{\sigma} \cdot \vec{n} = -\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^+}, \text{ and}\\
\boldsymbol{\sigma} \cdot \vec{n} = +\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^-}.
\end{gather}
%
The differentiation index is 2 because we must take the second time derivative of the prescribed slip equation to match the order of the time derivative in the elasticity equation,
%
\begin{gather}
\frac{\partial \vec{v}^+}{\partial t} - \frac{\partial \vec{v}^-}{\partial t} -
\frac{\partial^2 \vec{d}(\vec{x},t)}{\partial t^2} = \vec{0}.
\end{gather}
%
We generate the weak form in the usual way,
%
\begin{equation}
% Displacement-velocity
\int_{\Omega} {\vec{\psi}_\mathit{trial}^{u}} \cdot \frac{\partial \vec{u}}{\partial t} \, d\Omega = \int_{\Omega} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{v} \, d\Omega, \\
% Elasticity
\end{equation}
\begin{multline}
\int_{\Omega} {\vec{\psi}_\mathit{trial}^{v}} \cdot \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} \, d\Omega = \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\boldsymbol{\sigma}(\vec{u}) \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma \\
 + \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left(-\vec{\lambda}(\vec{x},t)\right) + {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left(+\vec{\lambda}(\vec{x},t)\right)\, d\Gamma, \\
\end{multline}
% Prescribed slip
\begin{equation}
\int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \left(\frac{\partial \vec{v}^+}{\partial t} - \frac{\partial \vec{v}^-}{\partial t} - \frac{\partial^2 \vec{d}(\vec{x},t)}{\partial t^2} \right) \, d\Gamma = 0.
\end{equation}
%
For compatibility with PETSc TS IMEX implementations, we need $\dot{\vec{s}}$ on the LHS for the explicit part (displacement-velocity and elasticity equations) and we need $\vec{\lambda}$ in the equation for the implicit part (prescribed slip equation).
We first focus on the explicit part and create a lumped LHS Jacobian matrix, $M$, so that we have

```{math}
:label: eqn:displacement:velocity:prescribed:slip:weak:form
% Displacement-velocity
\frac{\partial \vec{u}}{\partial t} = M_u^{-1} \int_{\Omega} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{v} \, d\Omega, \\
```
```{math}
:label: eqn:elasticity:prescribed:slip:dynamic:weak:form
% Elasticity
\begin{multline}
\frac{\partial \vec{v}}{\partial t} = M_v^{-1} \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\boldsymbol{\sigma}(\vec{u}) \, d\Omega + M_v^{-1} \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma \\
 + M_{v^+}^{-1} \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left(-\vec{\lambda}(\vec{x},t)\right) \, d\Gamma + M_{v^-}^{-1} \int_{\Gamma_{f}}{\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left(+\vec{\lambda}(\vec{x},t)\right) \, d\Gamma, \\
\end{multline}
```

\begin{gather}
% Mu
M_u = \mathit{Lump}\left( \int_\Omega {\psi_\mathit{trial}^{u}}_i \delta_{ij} {\psi_\mathit{basis}^{u}}_j \, d\Omega \right), \\
% Mv
M_v = \mathit{Lump}\left( \int_\Omega {\psi_\mathit{trial}^{v}}_i \rho(\vec{x}) \delta_{ij} {\psi_\mathit{basis}^{v}}_j \, d\Omega \right).
\end{gather}
%
Now, focusing on the implicit part we want to introduce $\vec{\lambda}$ into the prescribed slip equation.
We solve the elasticity equation for $\frac{\partial \vec{v}}{\partial t}$, create the weak form, and substitute into the prescribed slip equation.
Solving the elasticity equation for $\frac{\partial \vec{v}}{\partial t}$, we have
%
\begin{equation}
\frac{\partial \vec{v}}{\partial t} = \frac{1}{\rho(x)} \vec{f}(\vec{x},t) + \frac{1}{\rho(x)} \left(\nabla \cdot \boldsymbol{\sigma}(\vec{u}) \right),
\end{equation}
and the corresponding weak form is

```{math}
:label: eqn:prescribed:slip:DAE:weak:form
\int_{\Omega} {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{\partial \vec{v}}{\partial t} \, d\Omega = \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{1}{\rho(x)} \vec{f}(\vec{x},t) + {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \boldsymbol{\sigma}(\vec{u}) \right) \, d\Omega,
```

We apply the divergence theorem,
%
\begin{equation}
\int_{\Omega} \nabla \cdot \vec{F} \, d\Omega = \int_\Gamma \vec{F} \cdot \vec{n} \, d\Gamma,
\end{equation}
%
with $\vec{F} = {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \boldsymbol{\sigma}(\vec{u})\right)$ to get
%
\begin{equation}
\int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \boldsymbol{\sigma}(\vec{u}) \right) \, d\Omega = \int_\Gamma {\vec{\psi}_\mathit{trial}^{v}} \cdot \left( \frac{1}{\rho(x)} \boldsymbol{\sigma}(\vec{u}) \cdot \vec{n} \right) \, d\Gamma + \int_\Omega \nabla{\vec{\psi}_\mathit{trial}^{v}} : \left(-\frac{1}{\rho(\vec{x})} \boldsymbol{\sigma}(\vec{u}) \right) + {\vec{\psi}_\mathit{trial}^{v}} \cdot \left(-\frac{\nabla \rho(\vec{x})}{\rho^2} \cdot \boldsymbol{\sigma}(\vec{u}) \right) \, d\Omega.
\end{equation}
%
Restricting the trial function to $v^+$ and $v^-$ while recognizing that there is no overlap between the external Neumann boundary conditions $\Gamma_\tau$ and the fault interfaces $\Gamma_f$, yields
%
\begin{gather}
\int_\Omega {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \boldsymbol{\sigma}(\vec{u}) \right) \, d\Omega = \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left(-\frac{1}{\rho(x)} \vec{\lambda} \right) \, d\Gamma + \int_\Omega \nabla{\vec{\psi}_\mathit{trial}^{v^+}} : \left(-\frac{1}{\rho(\vec{x})} \boldsymbol{\sigma}(\vec{u}) \right) + {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \, \left(-\frac{\nabla \rho(\vec{x})}{\rho^2} \cdot \boldsymbol{\sigma}(\vec{u}) \right) \, d\Omega, \\
\int_\Omega {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \boldsymbol{\sigma}(\vec{u}) \right) \, d\Omega = \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left(+\frac{1}{\rho(x)} \vec{\lambda} \right) \, d\Gamma + \int_\Omega \nabla{\vec{\psi}_\mathit{trial}^{v^-}} : \left(-\frac{1}{\rho(\vec{x})} \boldsymbol{\sigma}(\vec{u}) \right) + {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \, \left(-\frac{\nabla \rho(\vec{x})}{\rho^2} \cdot \boldsymbol{\sigma}(\vec{u}) \right) \, d\Omega.
\end{gather}
%
Picking ${\vec{\psi}_\mathit{trial}^{v}}={\vec{\psi}_\mathit{trial}^{\lambda}}$ and substituting into equation {math:numref}`eqn:prescribed:slip:DAE:weak:form` gives
%
\begin{gather}
\int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \left( \frac{\partial \vec{v}^+}{\partial t} - \frac{\partial \vec{v}^-}{\partial t} - \frac{\partial^2 \vec{d}(\vec{x},t)}{\partial t^2} \right) \, d\Gamma = \\
\int_\Omega {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left( \frac{1}{\rho(\vec{x})} \vec{f}(\vec{x}, t) - \frac{\nabla \rho(\vec{x})}{\rho^2(\vec{x})} \cdot \boldsymbol{\sigma}(\vec{u}) \right) + \nabla {\vec{\psi}_\mathit{trial}^{v^+}} : \left(-\frac{1}{\rho(\vec{x})} \boldsymbol{\sigma}(\vec{u}) \right) \, d\Omega + \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left(-\frac{1}{\rho(\vec{x})} \vec{\lambda} \right) \, d\Gamma \\
- \int_\Omega {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left( \frac{1}{\rho(\vec{x})} \vec{f}(\vec{x}, t) -\frac{\nabla \rho(\vec{x})}{\rho^2(\vec{x})} \cdot \boldsymbol{\sigma}(\vec{u}) \right) + \nabla {\vec{\psi}_\mathit{trial}^{v^-}} : \left(-\frac{1}{\rho(\vec{x})} \boldsymbol{\sigma}(\vec{u}) \right) \, d\Omega + \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left(-\frac{1}{\rho(\vec{x})} \vec{\lambda} \right) \, d\Gamma \\
- \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \frac{\partial^2 \vec{d}(\vec{x}, t)}{\partial t^2} \, d\Gamma.
\end{gather}
%
We rewrite the integrals over the domain involving the degrees of freedom adjacent to the fault as integrals over the positive and negative sides of the fault.
These are implemented as integrals over the faces of cells adjacent to the fault; they involve quantities, such as density, that are defined only within the domain cells.
After collecting and rearranging terms, we have

```{math}
:label: eqn:elasticity:prescribed:slip:dynamic:DAE:weak:form
\begin{gather}
\int_{\Gamma_{f^+}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \frac{1}{\rho(\vec{x})} \left(\vec{\lambda} - \vec{f}(\vec{x},t) + \frac{\nabla\rho(\vec{x})}{\rho(\vec{x})} \cdot \boldsymbol{\sigma}(\vec{u}) \right) + \nabla {\vec{\psi}_\mathit{trial}^{\lambda}} : \left(+\frac{1}{\rho(\vec{x})} \boldsymbol{\sigma}(\vec{u}) \right) \, d\Gamma \\
+ \int_{\Gamma_{f^-}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \frac{1}{\rho(\vec{x})} \left(\vec{\lambda} + \vec{f}(\vec{x},t) - \frac{\nabla\rho(\vec{x})}{\rho(\vec{x})} \cdot \boldsymbol{\sigma}(\vec{u}) \right) + \nabla {\vec{\psi}_\mathit{trial}^{\lambda}} : \left(-\frac{1}{\rho(\vec{x})} \boldsymbol{\sigma}(\vec{u}) \right)  \, d\Gamma \\
+ \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \frac{\partial^2 \vec{d}(\vec{x}, t)}{\partial t^2} \, d\Gamma = 0.
\end{gather}
```

## Residual Pointwise Functions

Combining the explicit parts of the weak form in equations {math:numref}`eqn:displacement:velocity:prescribed:slip:weak:form` and {math:numref}`eqn:elasticity:prescribed:slip:dynamic:weak:form` with the implicit part of the weak form in equation {math:numref}`eqn:elasticity:prescribed:slip:dynamic:DAE:weak:form` and identifying $F(t,s,\dot{s})$ and $G(t,s)$, we have
%
\begin{gather}
% Fu
F^u(t,s,\dot{s}) = \frac{\partial \vec{u}}{\partial t} \\
% Fv
F^v(t,s,\dot{s}) = \frac{\partial \vec{v}}{\partial t} \\
% Fl
\end{gather}
\begin{multline}
F^\lambda(t,s,\dot{s}) =   \int_{\Gamma_{f^+}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot{\color{blue}\underbrace{\color{black}\frac{1}{\rho(\vec{x})} \left(  \vec{\lambda} - \vec{f}(\vec{x},t) + \frac{\nabla\rho(\vec{x})}{\rho(\vec{x})} \cdot \boldsymbol{\sigma}(\vec{u}) \right)}_{\color{blue}{f^\lambda_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{\lambda}} : {\color{blue}\underbrace{\color{black}\left(+\frac{1}{\rho(\vec{x})} \boldsymbol{\sigma}(\vec{u})\right)}_{\color{blue}{f^\lambda_1}}} \, d\Gamma \\
+ \int_{\Gamma_{f^-}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot {\color{blue} \underbrace{\color{black}\frac{1}{\rho(\vec{x})} \left(\vec{\lambda} + \vec{f}(\vec{x},t) - \frac{\nabla\rho(\vec{x})}{\rho(\vec{x})} \cdot \boldsymbol{\sigma}(\vec{u}) \right)}_{\color{blue}{f^\lambda_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{\lambda}} : {\color{blue}\underbrace{\color{black}\left(-\frac{1}{\rho(\vec{x})} \boldsymbol{\sigma}(\vec{u}) \right)}_{\color{blue}{f^\lambda_1}}}  \, d\Gamma \\
+ \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot {\color{blue}\underbrace{\color{black}\frac{\partial^2 \vec{d}(\vec{x}, t)}{\partial t^2}}_{\color{blue}{f^\lambda_0}}} \, d\Gamma \\
\end{multline}
% Gu
\begin{gather}
G^u(t,s) = \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot{\color{blue}\underbrace{\color{black}\vec{v}}_{\color{blue}{\vec{g}^u_0}}} \, d\Omega, \\
% Gv
G^v(t,s) =  \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot{\color{blue}\underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{g}^v_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{v}} : {\color{blue}\underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u})}_{\color{blue}{\boldsymbol{g^v_1}}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot {\color{blue}\underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{g}^v_0}}} \, d\Gamma, + \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot {\color{blue}\underbrace{\color{black}\left(-\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}{\vec{g}^v_0}}} + {\vec{\psi}_\mathit{trial}^{v^-}} \cdot {\color{blue}\underbrace{\color{black}\left(+\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}{\vec{g}^v_0}}} \, d\Gamma, \\
% Gl
G^l(t,s) = 0
\end{gather}
%
## Jacobian Pointwise Functions

For the explicit part we have pointwise functions for computing the lumped LHS Jacobian.
These are exactly the same pointwise functions as in the dynamic case without a fault,
%
\begin{equation}
\begin{aligned}
% J_F uu
J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{u}} = \int_\Omega {\psi_\mathit{trial}^{u}}_i{\color{blue}\underbrace{\color{black}s_\mathit{tshift} \delta_{ij}}_{\color{blue}{J^{uu}_{f0}}}} {\psi_\mathit{basis}^{u}}_j  \, d\Omega, \\
% J_F vv
J_F^{vv} &= \frac{\partial F^v}{\partial v} + s_\mathit{tshift} \frac{\partial F^v}{\partial \dot{v}} = \int_\Omega {\psi_\mathit{trial}^{v}}_i{\color{blue}  \underbrace{\color{black}\rho(\vec{x}) s_\mathit{tshift} \delta_{ij}}_{\color{blue}{J ^{vv}_{f0}}}} {\psi_\mathit{basis}^{v}}_j \, d\Omega
\end{aligned}
\end{equation}
%
For the implicit part, we have pointwise functions for the LHS Jacobians associated with the prescribed slip,
%
\begin{multline}
% J_F lu
J_F^{\lambda u} = \frac{\partial F^\lambda}{\partial u} + s_\mathit{tshift} \frac{\partial F^\lambda}{\partial \dot{u}} = \\
\int_{\Gamma_{f^+}} {\psi_\mathit{trial}^{\lambda}}_i {\color{blue}\underbrace{\color{black}\frac{\rho_{,j}(\vec{x})}{\rho^2(\vec{x})} C_{ikjl} {\psi_\mathit{basis}^{u}}_{j,l}}_{\color{blue}{J^{\lambda u}_{fX}}}} + {\psi_\mathit{trial}^{\lambda}}_{i,k}{\color{blue}  \underbrace{\color{black}+\frac{1}{\rho(\vec{x})} C_{ikjl} {\psi_\mathit{basis}^{u}}_{k,l}}_{\color{blue}{J^{\lambda u}_{f3}}}} \, d\Gamma \\
+\int_{\Gamma_{f^-}} {\psi_\mathit{trial}^{\lambda}}_i{\color{blue}\underbrace{\color{black}-\frac{\rho_{,j}(\vec{x})}{\rho^2(\vec{x})} C_{ikjl} {\psi_\mathit{basis}^{u}}_{j,l}}_{\color{blue}{J^{\lambda u}_{fX}}}} + {\psi_\mathit{trial}^{\lambda}}_{i,k} {\color{blue}\underbrace{\color{black}-\frac{1}{\rho(\vec{x})} C_{ikjl} {\psi_\mathit{basis}^{u}}_{k,l}}_{\color{blue}{J^{\lambda u}_{f3}}}} \, d\Gamma \\
\end{multline}
% J_F ll
\begin{equation}
J_F^{\lambda \lambda} = \frac{\partial F^\lambda}{\partial \lambda} + s_\mathit{tshift} \frac{\partial F^\lambda}{\partial \dot{\lambda}} = \int_{\Gamma_{f^+}} {\psi_\mathit{trial}^{\lambda}}_i {\color{blue}\underbrace{\color{black}\frac{1}{\rho(\vec{x})} \delta_{ij}}_{\color{blue}{J^{\lambda\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j \, d\Gamma + \int_{\Gamma_{f^-}} {\psi_\mathit{trial}^{\lambda}}_i {\color{blue}\underbrace{\color{black}\frac{1}{\rho(\vec{x})} \delta_{ij}}_{\color{blue}{J^{\lambda\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j \, d\Gamma
\end{equation}
