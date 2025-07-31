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
\frac{\partial \vec{v}}{\partial t} = M_v^{-1} \left(\int_\Omega \vec{\psi}_\mathit{trial}^{v} \cdot \vec{f}(\vec{x},t) + \nabla \vec{\psi}_\mathit{trial}^{v} : -\boldsymbol{\sigma}(\vec{u}) \, d\Omega  + \int_{\Gamma_\tau} \vec{\psi}_\mathit{trial}^{v} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma \right), \\
% Mu
M_u = \mathit{Lump}\left( \int_\Omega {\psi_\mathit{trial}^{u}}_i \delta_{ij} {\psi_\mathit{basis}^{u}}_j \, d\Omega \right), \\
% Mv
M_v = \mathit{Lump}\left( \int_\Omega {\psi_\mathit{trial}^{v}}_i \rho(\vec{x}) \delta_{ij} {\psi_\mathit{basis}^{v}}_j \, d\Omega \right).
\end{gather}

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

(sec-user-governing-eqns-absorbing-boundary)=
## Absorbing Boundary

The absorbing boundary is a special case of a Neumann boundary condition, in which the traction depends on the deformation. Consider a plane wave propagating at velocity $c$.
We can write the displacement field as
%
\begin{equation}
\vec{u}(\vec{x},t) = \vec{u^t}(t-\frac{\vec{x}}{c}),
\end{equation}
%
where $\vec{x}$ is position, $t$ is time, and $\vec{u^t}$ is the shape of the propagating wave.
For an absorbing boundary we want the traction on the boundary to be equal to the traction associated with the wave propagating out of the domain.
Starting with the expression for the traction on a boundary, $T_{i}=\sigma_{ij}n_{j},$ and using the local coordinate system for the boundary $s_{h}s_{v}n,$ where $\vec{n}$ is the direction normal to the boundary, $\vec{s}_h$ is the horizontal direction tangent to the boundary, and $\vec{s}_v$ is the vertical direction tangent to the boundary, the tractions on the boundary are
%
\begin{align}
\tau_{s_h} &= \sigma_{s_{h}n}\\
\tau_{s_v} &= \sigma_{s_{v}n}\\
\tau_{n} &=\sigma_{nn}.
\end{align}

In the case of a horizontal boundary, we define a reference directions in order to assign unique tangential directions.
For a linear elastic isotropic material, $\sigma_{ij}=\lambda\epsilon_{kk}\delta_{ij}+2\mu\epsilon_{ij},$ and we can write the tractions as
%
\begin{align}
\tau_{s_{h}} &= 2 \mu \epsilon_{s_{h}n}\\
\tau_{s_{v}} &= 2 \mu \epsilon_{s_{v}n}\\
\tau_{n} &= (\lambda+2\mu) \epsilon_{nn} + \lambda (\epsilon_{s_{h}s_{h}} + \epsilon_{s_{v}s_{v}}).
\end{align}
%
For infinitesimal strains, $\epsilon_{ij}=\frac{1}{2}(u_{i,j}+u_{j,i})$ and we have
%
\begin{align}
\epsilon_{s_{h}n} &= \frac{1}{2} (u_{s_{h},n} + u_{n,s_{h}})\\
\epsilon_{s_{v}n} &= \frac{1}{2} (u_{s_{v},n} + u_{n,s_{v}})\\
\epsilon_{nn} &= u_{n,n}.
\end{align}
%
For our propagating plane wave, we recognize that
%
\begin{equation}
\frac{\partial\vec{u^t}(t-\frac{\vec{x}}{c})}{\partial x_{i}} = -\frac{1}{c} \frac{\partial\vec{u^{t}}(t-\frac{\vec{x}}{c})}{\partial t} = -\frac{1}{c} \vec{v}^t(t-\frac{\vec{x}}{c}),
\end{equation}
%
so that our expressions for the tractions become
%
\begin{gather}
\tau_{s_{h}} = -\frac{\mu}{c} \left(v_{s_{h}}^{t}(t-\frac{\vec{x}}{c})+v_n^t(t-\frac{\vec{x}}{c})\right),\\
\tau_{s_{v}} = -\frac{\mu}{c} \left(v_{s_{v}}^{t}(t-\frac{\vec{x}}{c})+v_{n}^{t}(t-\frac{\vec{x}}{c})\right).
\end{gather}

For the normal traction consider a dilatational wave propagating normal to the boundary at speed $v_p$; in this case $u_{s_{h}}=u_{s_{v}}=0$ and $c=v_{p}$.
For the shear tractions, consider a shear wave propagating normal to the boundary at speed $v_s$; we can decompose this into one case where $u_{n}=u_{s_{v}}=0$ and another case where $u_{n}=u_{s_{h}}=0$, with $c=v_{s}$ in both cases.
We also recognize that $\mu=\rho v_{s}^{2}$ and $\lambda+2\mu=\rho v_{p}^{2}$.
This leads to the following expressions for the tractions:
%
\begin{gather}
\tau_{s_{h}}=-\rho v_{s} v_{s_{h}}^{t}(t-\frac{\vec{x}}{c})\\
\tau_{s_{v}}=-\rho v_{s} v_{v}^{t}(t-\frac{\vec{x}}{c})\\
\tau_{n}=-\rho v_{p} v_{n}^{t}(t-\frac{\vec{x}}{c})
\end{gather}

Substituting the tractions into the weak form for the Neumann boundary, we have
%
\begin{equation}
\int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma = \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \left( -\rho c_i v_i(t) \right) \, d\Gamma,
\end{equation}
%
where $c_i$ equals $v_p$ for the normal traction and $v_s$ for the shear tractions.

### Residual Pointwise Function

\begin{equation}
G^v(t,s) =  \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot{\color{blue}\underbrace{\color{black} \left( -\rho c_i v_i(t) \right) }_{\color{blue}{\vec{g}^v_0}}} \, d\Gamma,
\end{equation}
