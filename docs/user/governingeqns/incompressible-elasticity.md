(sec-user-governing-eqns-incompressible-elasticity)=
# Incompressible Isotropic Elasticity with Infinitesimal Strain (Bathe)

In this section we apply a similar approach to the one we use for the elasticity equation to the case of an incompressible material.
We only consider the quasistatic case (neglect inertia) without faults.
As the bulk modulus ($K$) approaches infinity, the volumetric strain ($\mathop{\mathrm{Tr}}(\epsilon)$) approaches zero and the pressure remains finite, $p = -K \mathop{\mathrm{Tr}}(\epsilon)$.
We consider pressure $p$ as an independent variable and decompose the stress into the pressure and deviatoric components.
As a result, we write the stress tensor in terms of both the displacement and pressure fields,
%
\begin{equation}
\boldsymbol{\sigma}(\vec{u},p) = \boldsymbol{\sigma}^\mathit{dev}(\vec{u}) - p\boldsymbol{I}.
\end{equation}
%
The strong form is
%
\begin{gather}
  % Solution
  \vec{s}^T = \left( \vec{u} \quad \ p \right)^T, \\
  % Elasticity
  \vec{f}(t) + \boldsymbol{\nabla} \cdot \left(\boldsymbol{\sigma}^\mathit{dev}(\vec{u}) - p\boldsymbol{I}\right) = \vec{0} \text{ in }\Omega, \\
  % Pressure
  \vec{\nabla} \cdot \vec{u} + \frac{p}{K} = 0 \text{ in }\Omega, \\
  % Neumann
  \boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau} \text{ on }\Gamma_\tau, \\
  % Dirichlet
  \vec{u} = \vec{u}_0 \text{ on }\Gamma_u, \\
  p = p_0 \text{ on }\Gamma_p.
\end{gather}
%
We place all terms for the elasticity and pressure equations on the left-hand-side, consistent with PETSc TS implicit time stepping.

```{table} Mathematical notation for incompressible elasticity with infinitesimal strain
:name: tab:notation:incompressible:elasticity

| Category                       |               Symbol               | Description                                               |
| :----------------------------- | :--------------------------------: | :-------------------------------------------------------- |
| Unknowns                       |             $\vec{u}$              | Displacement field                                        |
|                                |                $p$                 | Pressure field ($p>0$ corresponds to negative mean stress |
| Derived quantities             |       $\boldsymbol{\sigma}$        | Cauchy stress tensor                                      |
|                                | $\boldsymbol{\sigma}^\mathit{dev}$ | Cauchy deviatoric stress tensor                           |
|                                |      $\boldsymbol{\epsilon}$       | Cauchy strain tensor                                      |
| Common constitutive parameters |               $\rho$               | Density                                                   |
|                                |               $\mu$                | Shear modulus                                             |
|                                |                $K$                 | Bulk modulus                                              |
| Source terms                   |             $\vec{f}$              | Body force per unit volume, for example $\rho \vec{g}$    |
```

Using trial functions ${\vec{\psi}_\mathit{trial}^{u}}$ and ${\psi_\mathit{trial}^{p}}$ and incorporating the Neumann boundary conditions, we write the weak form as
%
\begin{gather}
% Displacement
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{f}(t) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : \left(-\boldsymbol{\sigma}^\mathit{dev}(\vec{u}) + p\boldsymbol{I}
\right)\, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{\tau}(t) \, d\Gamma, = 0 \\
% Pressure
\int_\Omega {\psi_\mathit{trial}^{p}} \cdot \left(\vec{\nabla} \cdot \vec{u} + \frac{p}{K} \right) \, d\Omega = 0.
\end{gather}
%
## Residual Pointwise Functions

Identifying $F(t,s,\dot{s})$, we have
```{math}
:label: eqn:incompressible:elasticity:displacement
\begin{gather}
F^u(t,s,\dot{s}) = \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot{\color{blue}
\underbrace{\color{black}\vec{f}(t)}_{\color{blue}{f_0^u}}} + \nabla {\vec{\psi}_\mathit{trial}^{u}} :{\color{blue}
\underbrace{\color{black}\left(-\boldsymbol{\sigma}^\mathit{dev}(\vec{u}) + p\boldsymbol{I}\right)}_{\color{blue}{f_1^u}}}  \, d\Omega
+ \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}\underbrace{\color{black}\vec{\tau}(t)}_{\color{blue}{f_0^u}}} \, d\Gamma, \\
%
\end{gather}
```
```{math}
:label: eqn:incompressible:elasticity:pressure
\begin{gather}
F^p(t,s,\dot{s}) = \int_\Omega {\psi_\mathit{trial}^{p}} \cdot {\color{blue}\underbrace{\color{black}\left(\vec{\nabla} \cdot \vec{u} +
\frac{p}{K} \right)}_{\color{blue}{f_0^p}}} \, d\Omega.
\end{gather}
```
## Jacobians Pointwise Functions

With two fields we have four Jacobian pointwise functions for the LHS:
%
\begin{align}
% JF uu
J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{u}} =
           \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial u}(-\boldsymbol{\sigma}^\mathit{dev}) \, d\Omega
            = \int_\Omega {\psi_\mathit{trial}^{u}}_{i,k} \, {\color{blue}
\underbrace{\color{black}\left(-C^\mathit{dev}_{ikjl}\right)}_{\color{blue}{J_{f3}^{uu}}}}  \, {\psi_\mathit{basis}^{u}}_{j,l}\, d\Omega \\
% JF up
J_F^{up} &= \frac{\partial F^u}{\partial p} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{p}} =
           \int_\Omega \nabla{\vec{\psi}_\mathit{trial}^{u}} : \boldsymbol{I} {\psi_\mathit{basis}^{p}} \,  d\Omega
           = \int_\Omega {\psi_\mathit{trial}^{u}}_{i,k} {\color{blue}  \underbrace{\color{black}\delta_{ik}}_{\color{blue}{J_{f2}^{up}}}} \, {\psi_\mathit{basis}^{p}} \, d\Omega \\
% JF pu
J_F^{pu} &= \frac{\partial F^p}{\partial u} + s_\mathit{tshift} \frac{\partial F^p}{\partial \dot{u}} =
           \int_\Omega {\psi_\mathit{trial}^{p}} \left(\vec{\nabla}  \cdot {\vec{\psi}_\mathit{basis}^{u}}\right) \, d\Omega
           = \int_\Omega {\psi_\mathit{trial}^{p}} {\color{blue}\underbrace{\color{black}\delta_{jl}}_{\color{blue}{J_{f1}^{pu}}}} {\psi_\mathit{basis}^{u}}_{j,l} \, d\Omega\\
% JF pp
J_F^{pp} &= \frac{\partial F^p}{\partial p}  + s_\mathit{tshift} \frac{\partial F^p}{\partial \dot{p}} =
           \int_\Omega {\psi_\mathit{trial}^{p}}{\color{blue}\underbrace{\color{black}\frac{1} {K}}_{\color{blue}{J_{f0}^{pp}}}} {\psi_\mathit{basis}^{p}} \, d\Omega
\end{align}
%
For isotropic, linear incompressible elasticity, the deviatoric elastic constants are:
%
\begin{align}
  C_{1111} &= C_{2222} = C_{3333} = +\frac{4}{3} \mu \\
  C_{1122} &= C_{1133} = C_{2233} = -\frac{2}{3} \mu \\
  C_{1212} &= C_{1313} = C_{2323} = \mu
\end{align}
