# Elasticity with Infinitesimal Strain and Prescribed Slip on Faults

For each fault, which is an internal interface, we add a boundary condition to the elasticity equation prescribing the jump in the displacement field across the fault,
%
```{math}
:label: eqn:bc:prescribed:slip
\begin{gather}
  \vec{u}^+ - \vec{u}^- - \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f,
\end{gather}
```
%
where $\vec{u}^+$ is the displacement vector on the "positive" side of the fault, $\vec{u}^-$ is the displacement vector on the "negative" side of the fault, $\vec{d}$ is the slip vector on the fault, and $\vec{n}$ is the fault normal which points from the negative side of the fault to the positive side of the fault.
We enforce the jump in displacements across the fault using a Lagrange multiplier corresponding to equal and opposite tractions on the two sides of the fault.

We apply conservation of momemtum,
\begin{equation}
  \int_\Omega \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} \, d\Omega = \int_\Omega \vec{f}(\vec{x},t) \, d\Omega + \int_\Gamma \vec{\tau}(\vec{x},t) \, d\Gamma,
\end{equation}
to a fault interface $\Omega_f$ with boundaries $\Gamma_{f^+}$ and $\Gamma_{f^-}$.
For a fault interface, the body force is zero, $\vec{f}(\vec{x},t) = \vec{0}$.
The tractions on the positive and negative fault faces are
\begin{gather}
  \tau^+(\vec{x},t) = \boldsymbol{\sigma}^+ \cdot \vec{n} + \vec{\lambda} \\
  \tau^-(\vec{x},t) = \boldsymbol{\sigma}^- \cdot \vec{n} - \vec{\lambda},
\end{gather}

where $\vec{\lambda}$ is the Lagrange multiplier that corresponds to the fault traction generating the prescribed slip and $\boldsymbol{\sigma}^+$ and $\boldsymbol{\sigma}^-$ are the stresses in the domain at the positive and negative sides of the fault.
Thus, for a fault interface, we have
\begin{equation}
  \int_{\Omega_f} \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} \, d\Omega = \int_{\Gamma_{f^+}} \boldsymbol{\sigma} \cdot \vec{n} + \vec{\lambda} \, d\Gamma + \int_{\Gamma_{f^-}} \boldsymbol{\sigma} \cdot \vec{n} - \vec{\lambda} \, d\Gamma.
\end{equation}

```{table} Mathematical notation for elasticity equation with infinitesimal strain and prescribed slip on faults.
:name: tab:notation:elasticity:prescribed:slip
| Category | Symbol | Description |
|:-------------|:----------:| :-------------------------|
| Unknowns |    $\vec{u}$    | Displacement field       |
|          |    $\vec{v}$    | Velocity field           |
|          | $\vec{\lambda}$ | Lagrange multiplier field |
| Derived quantities |  $\boldsymbol{\sigma}$  | Cauchy stress tensor |
|                    | $\boldsymbol{\epsilon}$ | Cauchy strain tensor |
| Common constitutive parameters | $\rho$  | Density        |
|                                | $\mu$   | Shear modulus  |
|                                | $K$     | Bulk modulus   |
| Source terms  | $\vec{f}$ | Body force per unit volume, for example $\rho \vec{g}$ |
|               | $\vec{d}$ | Slip vector field on the fault corresponding to a jump in the displacement field across the fault |
```

:::{toctree}
quasistatic.md
dynamic.md
:::
