# Nondimensionalization of Elasticity Equation

Starting with the elasticity equation,
%
\begin{gather}
\rho(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}=\vec{0}\text{ in }\Omega,\\
\boldsymbol{\sigma}(\vec{u})\cdot\vec{n}=\vec{\tau}(\vec{x},t)\text{ on }\Gamma_{\tau}\text{,}\\
\vec{u}^{+}-\vec{u}^{-}=\vec{d}\text{ on }\Gamma_{f}.
\end{gather}
%
we define nondimensional values:
%
\begin{align}
\vec{x}^* &= \frac{\vec{x}}{x_o}, \\
\vec{u}^* &= \frac{\vec{u}}{u_o}, \\
\rho^* &= \frac{\rho}{\rho_o}, \\
\vec{f}^* &= \frac{\vec{f}}{f_o}, \\
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
\rho_o  \rho^*(\vec{x}^*) \frac{u_o}{t_o^2} \frac{\partial^{2}\vec{u}^*}{\partial {t^*}^{2}} - f_o \vec{f}^*(\vec{x}^*,t^*) - \frac{\sigma_o}{x_o} \boldsymbol{\nabla}^*\cdot\boldsymbol{\sigma}^* = \vec{0}\text{ in }\Omega,\\
\sigma_o \boldsymbol{\sigma}^*(\vec{u^*}) \cdot \vec{n} = \sigma_o \vec{\tau}^*(\vec{x}^*,t^*)\text{ on }\Gamma_{\tau}\text{,}\\
u_o \left(\vec{u}^{*^{+}}-\vec{u}^{*^{-}}\right) = u_o \vec{d}^* \text{ on }\Gamma_{f}.
\end{gather}

For the second two equations, the nondimensional scales in for the terms in each equation are consistent, so we will limit our discussion to the first equation.
Grouping terms and multiplying by $\frac{x_o}{\sigma_o}$, we have
%
\begin{equation}
\left( \rho_o \frac{u_o}{t_o^2}\frac{x_o}{\sigma_o}\right) \rho^*(\vec{x}^*) \frac{\partial^{2}\vec{u}^*}{\partial {t^*}^{2}} - \left(f_o \frac{x_o}{\sigma_o}\right) \vec{f}^*(\vec{x}^*,t^*) -  \boldsymbol{\nabla}^*\cdot\boldsymbol{\sigma}^* = \vec{0}\text{ in }\Omega.
\end{equation}
%
All terms should be nondimensional, which implies
%
\begin{align}
f_o &= \frac{\sigma_o}{x_o}, \\
\rho_o &= \frac{t_o^2 \sigma_o}{u_o x_o}.
\end{align}

We want to determine the stress scale, $\sigma_o$.
Considering isotropic, linear elasticity we have
%
\begin{equation}
\boldsymbol{\sigma} = \boldsymbol{C} : \boldsymbol{\epsilon} = \boldsymbol{C} : \frac{1}{2}\left(\boldsymbol{\nabla} + \boldsymbol{\nabla}^T \right) \vec{u}.
\end{equation}
%
Substituting in our nondimensional values yields
%
\begin{equation}
\sigma_o \boldsymbol{\sigma}^* = \mu_o \boldsymbol{C}^* : \frac{u_o}{x_o} \frac{1}{2}\left(\boldsymbol{\nabla}^* + \boldsymbol{\nabla}^{*^T}\right) \vec{u}^*.
\end{equation}
%
We recognize that for the equation to be nondimensional,
\begin{equation}
\sigma_o = \mu_o \frac{u_o}{x_o}.
\end{equation}

Returning to the expression for $\rho_o$ and substituting in the expression for $\sigma_o$, we have
\begin{align}
\rho_o &= \frac{t_o^2 \sigma_o}{u_o x_o}, \\
\rho_o &= \mu_o \frac{t_o^2}{x_o^2}.
\end{align}

## Inertia

The scale of the inertial term is $\Pi_\mathit{inertia} = \frac{\rho_o x_o^2}{\mu_o t_o^2}$.
If the time scale equals the time it takes the shear wave ($v_o^2 = \frac{\mu_o}{\rho_o}$) to propagate over the length scale, then
\begin{align}
t_o &= \frac{x_o}{v_o}, \\
t_o^2 &= \frac{x_o^2}{v_o^2}, \\
t_o^2 &= \frac{x_o^2 \rho_o}{\mu_o}.
\end{align}
Substituting into the expression for $\Pi_\mathit{inertia}$, we have
\begin{align}
\Pi_\mathit{inertia} &= \frac{\rho_o x_o^2}{\mu_o t_o^2}, \\
\Pi_\mathit{inertia} &= \frac{\rho_o x_o^2}{\mu_o} \frac{\mu_o}{\rho_o x_o^2}, \\
\Pi_\mathit{inertia} &= 1.
\end{align}
In this case, the inertial term is important, and we have the dynamic case with propagating seismic waves.

If the time scale is much larger than the time it takes the shear wave to propagate over the length scale, then
\begin{align}
t_o &\gg \frac{x_o}{v_o}, \\
t_o^2 &\gg \frac{x_o^2 \rho_o}{\mu_o}.
\end{align}
Substituting into the expression for $\Pi_\mathit{inertia}$, we have
\begin{align}
\Pi_\mathit{inertia} &= \frac{\rho_o x_o^2}{\mu_o t_o^2}, \\
\Pi_\mathit{inertia} &\ll \frac{\rho_o x_o^2}{\mu_o} \frac{\mu_o}{\rho_o x_o^2}, \\
\Pi_\mathit{inertia} &\ll 1.
\end{align}
In this case, the inertial term is negligible, and we have quasistatic elasticity.
:::

:::{note}
In PyLith v4 and earlier, we used a displacement scale equal to the length scale.
Using separate length and displacement scales facilitates accurately solving for displacements many orders of magnitude smaller than the length scale.
:::

```{table} Nondimensionalization of elasticity
:name: tab:elasticity:scales
| **Quantity** | **Scale**          |
| :---------: | :------------------ |
| $x_o$       | Length scale        |
| $u_o$       | Displacement scale  |
| $\mu_o$     | Rigidity scale      |
| $t_o$       | Time scale          |
| $\sigma_o = \mu_o \frac{u_o}{x_o}$  | Stress scale        |
| $f_o = \frac{\sigma_o}{x_o}$    | Density scale       |
| $\rho_o = \mu_o \frac{t_o^2}{x_o^2}$    | Density scale       |
```
