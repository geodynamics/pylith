# Nondimensionalization

Starting with the equations governing poroelasticity,
%
\begin{gather}
% Elasticity
  \rho_s\frac{\partial^2 \vec{u}}{\partial t^2} - \vec{f}(t) - \nabla \cdot \boldsymbol{\sigma}(\vec{u},p) = \vec{0} \text{ in } \Omega, \\
% Fluid mass balance
  \frac{\partial \zeta(\vec{u},p)}{\partial t} + \nabla \cdot \vec{q}(p) - \gamma(\vec{x},t) = 0 \text{ in } \Omega, \\
% Darcy flow
  \vec{q}(p) = -\frac{\boldsymbol{k}}{\mu_{f}}(\nabla p - \vec{f}_f),
\end{gather}
%
we define nondimensional values:
%
\begin{align}
\vec{x}^* &= \frac{\vec{x}}{x_o}, \\
\vec{u}^* &= \frac{\vec{u}}{u_o}, \\
p^* &= \frac{p}{p_o}, \\
\rho_s^* &= \frac{\rho_s}{\rho_o}, \\
\rho_f^* &= \frac{\rho_f}{\rho_o}, \\
\vec{f}^* &= \frac{\vec{f}}{f_o}, \\
\zeta^* &= \frac{\zeta}{\zeta_o}, \\
\vec{q}* &= \frac{\vec{q}}{q_o}, \\
\gamma* &= \frac{\gamma}{\gamma_o}, \\
\boldsymbol{\sigma}^* &= \frac{\boldsymbol{\sigma}}{\sigma_o}, \\
\vec{\tau}^* &= \frac{\vec{\tau}}{\sigma_o}.
\end{align}
%
We also recognize that
%
\begin{equation}
\boldsymbol{\nabla}^* = x_o \boldsymbol{\nabla}.
\end{equation}

Substituting into the equations, we have
%
\begin{gather}
% Elasticity
\rho_o  \rho^*(\vec{x}^*) \frac{u_o}{t_o^2} \frac{\partial^{2}\vec{u}^*}{\partial {t^*}^{2}} - f_o \vec{f}^*(\vec{x}^*,t^*) - \frac{\sigma_o}{x_o} \boldsymbol{\nabla}^*\cdot\boldsymbol{\sigma}^* = \vec{0}\text{ in }\Omega,\\
% Fluid mass balance
\frac{\zeta_o}{t_o} \frac{\partial \zeta^*(\vec{u},p)}{\partial t^*} + \frac{q_o}{x_o} \left( \nabla^* \cdot \vec{q^*}(p)\right) - \gamma_o \gamma^*(\vec{x},t) = 0 \text{ in } \Omega, \\
% Darcy flow
q_o \vec{q^*}(p) = -\frac{\boldsymbol{k}}{\mu_{f}}(\frac{p_o}{x_o} \nabla^* p^* - f_o \vec{f^*}_f).
\end{gather}

Grouping terms we have
%
\begin{gather}
% Elasticity
\left( \rho_o \frac{u_o}{t_o^2}\frac{x_o}{\sigma_o}\right) \rho^*(\vec{x}^*) \frac{\partial^{2}\vec{u}^*}{\partial {t^*}^{2}} - \left(f_o \frac{x_o}{\sigma_o}\right) \vec{f}^*(\vec{x}^*,t^*) -  \boldsymbol{\nabla}^*\cdot\boldsymbol{\sigma}^* = \vec{0}\text{ in }\Omega, \\
% Fluid mass balance
\left(\frac{\zeta_o}{t_o}\right) \frac{\partial \zeta^*(\vec{u},p)}{\partial t^*} + \left(\frac{q_o}{x_o}\right) \left( \nabla^* \cdot \vec{q^*}(p)\right) - \gamma_o \gamma^*(\vec{x},t) = 0 \text{ in } \Omega, \\
% Darcy flow
q_o \vec{q^*}(p) = -\left(\frac{p_o}{x_o} \frac{\boldsymbol{k}}{\mu_{f}}\right) ( \nabla^* p^* - f_o \frac{x_o}{p_o} \vec{f^*}_f).
\end{gather}

We also have
\begin{equation}
  \zeta(\vec{u},p) = \alpha (\nabla \cdot \vec{u}) + \frac{p}{M},
\end{equation}
%
and substituting in the nondimensionalized values leads to
%
\begin{equation}
  \zeta_o \zeta^*(\vec{u},p) = \alpha \frac{u_o}{x_o} (\nabla^* \cdot \vec{u}^*) + \frac{p_o}{\mu_o} \frac{p^*}{M^*},
\end{equation}
where $M^* = \frac{M}{\mu_o}$.
%
Grouping terms and recognizing that $p_o = \mu_o \frac{u_o}{x_o}$, we have
\begin{equation}
  \zeta_o \zeta^*(\vec{u},p) = \left(\alpha \frac{u_o}{x_o}\right) (\nabla^* \cdot \vec{u}^*) + \frac{u_o}{x_o} \frac{p^*}{M^*}.
\end{equation}
All terms should be nondimensional, which implies
%
\begin{align}
\zeta_o &= \frac{u_o}{x_o}.
\end{align}

Returning to the fluid mass balance, we have
%
\begin{equation}
\left(\frac{u_o}{t_o x_o}\right) \frac{\partial \zeta^*(\vec{u},p)}{\partial t^*} + \left(\frac{q_o}{x_o}\right) \left( \nabla^* \cdot \vec{q^*}(p)\right) - \gamma_o \gamma^*(\vec{x},t) = 0 \text{ in } \Omega, \\
\end{equation}
%
All terms should be nondimensional, which implies
%
\begin{align}
q_o &= \frac{u_o}{t_o}, \\
\gamma_o &= \frac{1}{t_o} \frac{u_o}{x_o}.
\end{align}

Regrouping terms leads to
%
\begin{equation}
\left(\frac{u_o}{t_o x_o}\right) \frac{\partial \zeta^*(\vec{u},p)}{\partial t^*} + \left( \nabla^* \cdot \vec{q^*}(p)\right) - \gamma^*(\vec{x},t) = 0 \text{ in } \Omega.
\end{equation}

Returning to the equation for Darcy flow, we have
%
\begin{equation}
\frac{u_o}{t_o} \vec{q^*}(p) = -\left(\frac{\sigma_o}{x_o} \frac{\boldsymbol{k}}{\mu_{f}}\right) ( \nabla^* p^* - \vec{f^*}_f),
\end{equation}
where we have made use of $q_o = \frac{u_o}{t_o}$, $\sigma_o = p_o$ and $f_o = \frac{\sigma_o}{x_p}$.
Rearranging we have,
\begin{align}
\vec{q^*}(p) &= -\left(\frac{\sigma_o t_o}{u_o x_o} \frac{\boldsymbol{k}}{\mu_{f}}\right) ( \nabla^* p^* - \vec{f^*}_f), \\
\vec{q^*}(p) &= -\left(\frac{\mu_o t_o}{x_o^2} \frac{\boldsymbol{k}}{\mu_{f}}\right) ( \nabla^* p^* - \vec{f^*}_f).
\end{align}
%
This implies
%
\begin{align}
\boldsymbol{k^*} = \frac{\boldsymbol{k}}{x_o^2}, \\
\mu_f^* = \frac{\mu_f}{\mu_o t_o}.
\end{align}
%
We also find that the time scale is given by
\begin{equation}
t_o = \frac{\mu_{f_o}}{\boldsymbol{k}_o} \frac{x_o^2}{\mu_o}.
\end{equation}

:::{note}
In PyLith v4 and earlier, we used a displacement scale equal to the length scale.
Using separate length and displacement scales facilitates accurately solving for displacements many orders of magnitude smaller than the length scale.
It also separates the stress scale used to nondimensionalize stress, tractions, and fluid pressure from the rigidity scale.
:::

```{table} Nondimensionalization of poroelasticity/
:name: tab:poroelasticity:scales
|                         **Quantity**                          | **Scale**          |
| :-----------------------------------------------------------: | :----------------- |
|                             $x_o$                             | Length scale       |
|                             $u_o$                             | Displacement scale |
|                            $\mu_o$                            | Rigidity scale     |
| $t_o= \frac{\mu_{f_o}}{\boldsymbol{k}_o} \frac{x_o^2}{\mu_o}$ | Time scale         |
|           $\sigma_o = p_o = \mu_o \frac{u_o}{x_o}$            | Stress scale       |
|                 $f_o = \frac{\sigma_o}{x_o}$                  | Body force scale   |
|             $\rho_o = \mu_o \frac{t_o^2}{x_o^2}$              | Density scale      |
|                   $q_o = \frac{u_o}{t_o} $                    | Flow scale         |
|          $\gamma_o = \frac{1}{t_o} \frac{u_o}{x_o}$           | Source scale       |
|              $\mu_f^* = \frac{\mu_f}{\mu_o t_o}$              | Viscosity scale    |
|       $\boldsymbol{k}_o = \frac{\boldsymbol{k}}{x_o^2}$       | Permeability scale |
```
