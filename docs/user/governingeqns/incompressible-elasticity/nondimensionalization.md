# Nondimensionalization

Starting with the incompressible elasticity equation,
%
\begin{gather}
% Elasticity
\rho(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}=\vec{0}\text{ in }\Omega, \\
% Pressure
\vec{\nabla} \cdot \vec{u} + \frac{p}{K} = 0 \text{ in }\Omega, \text{ where}\\
% Stress
\boldsymbol{\sigma} = \boldsymbol{\sigma}^{\mathit{dev}} - p \boldsymbol{I}.
\end{gather}
%
we define nondimensional values:
%
\begin{align}
\vec{x}^* &= \frac{\vec{x}}{x_o}, \\
\vec{u}^* &= \frac{\vec{u}}{u_o}, \\
p^* &= \frac{p}{p_o}, \\
\rho^* &= \frac{\rho}{\rho_o}, \\
\vec{f}^* &= \frac{\vec{f}}{f_o}, \\
\boldsymbol{\sigma}^* &= \frac{\boldsymbol{\sigma}}{\sigma_o}, \\
\vec{\tau}^* &= \frac{\vec{\tau}}{\sigma_o}.
\end{align}
%
As in the case of the elasticity equation, we also recognize that
%
\begin{equation}
\boldsymbol{\nabla}^* = x_o \boldsymbol{\nabla}.
\end{equation}

Substituting into the equations, we have
%
\begin{gather}
% Elasticity
\rho_o  \rho^*(\vec{x}^*) \frac{u_o}{t_o^2} \frac{\partial^{2}\vec{u}^*}{\partial {t^*}^{2}} - f_o \vec{f}^*(\vec{x}^*,t^*) - \frac{\sigma_o}{x_o} \boldsymbol{\nabla}^*\cdot\boldsymbol{\sigma}^* = \vec{0}\text{ in }\Omega,\\
% Pressure
\frac{u_o}{x_o}\left( \vec{\nabla}^* \cdot \vec{u}^*\right) + \frac{p_o}{\mu_o} \frac{p^*}{K^*} = 0 \text{ in }\Omega.
\end{gather}

Grouping terms we have
%
\begin{gather}
% Elasticity
\left( \rho_o \frac{u_o}{t_o^2}\frac{x_o}{\sigma_o}\right) \rho^*(\vec{x}^*) \frac{\partial^{2}\vec{u}^*}{\partial {t^*}^{2}} - \left(f_o \frac{x_o}{\sigma_o}\right) \vec{f}^*(\vec{x}^*,t^*) -  \boldsymbol{\nabla}^*\cdot\boldsymbol{\sigma}^* = \vec{0}\text{ in }\Omega, \\
% Pressure
\vec{\nabla}^* \cdot \vec{u}^* + \frac{p_o u_o}{\mu_o x_o} \frac{p^*}{K^*} = 0 \text{ in }\Omega.
\end{gather}
%
All terms should be nondimensional, which implies
%
\begin{align}
f_o &= \frac{\sigma_o}{x_o}, \\
\rho_o &= \frac{t_o^2 \sigma_o}{u_o x_o}, \\
p_o &= \mu_o \frac{x_o}{u_o}.
\end{align}

:::{note}
In PyLith v4 and earlier, we used a displacement scale equal to the length scale.
Using separate length and displacement scales facilitates accurately solving for displacements many orders of magnitude smaller than the length scale.
It also separates the stress scale used to nondimensionalize stress, tractions, and pressure from the rigidity scale.
:::

```{table} Nondimensionalization of incompressible elasticity.
:name: tab:incompressible:elasticity:scales
| **Quantity** | **Scale**          |
| :---------: | :------------------ |
| $x_o$       | Length scale        |
| $u_o$       | Displacement scale  |
| $\mu_o$     | Rigidity scale      |
| $t_o$       | Time scale          |
| $\sigma_o = p_o = \mu_o \frac{u_o}{x_o}$  | Stress scale        |
| $f_o = \frac{\sigma_o}{x_o}$    | Body force scale       |
| $\rho_o = \mu_o \frac{t_o^2}{x_o^2}$    | Density scale       |
```
