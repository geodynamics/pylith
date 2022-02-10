# Elasticity with Infinitesimal Strain and Prescribed Slip on Faults

For each fault, which is an internal interface, we add a boundary condition to the elasticity equation prescribing the jump in the displacement field across the fault,

```{math}
:label: eqn:bc:prescribed_slip
\begin{gather}
\vec{u}^+(\vec{x},t) - \vec{u}^-(\vec{x},t) - \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f,
\end{gather}
```
where $\vec{u}^+$ is the displacement vector on the "positive" side of the fault, $\vec{u}^-$ is the displacement vector on the "negative" side of the fault, $\vec{d}$ is the slip vector on the fault, and $\vec{n}$ is the fault normal which points from the negative side of the fault to the positive side of the fault.
Note that as an alternative to prescribing the jump in displacement across the fault, we can also prescribe the jump in velocity or acceleration across the fault in terms of slip rate or slip acceleration, respectively.

Using a domain decomposition approach for constraining the fault slip yields additional Neumann-like boundary conditions on the fault surface,
%
\begin{gather}
\boldsymbol{\sigma} \cdot \vec{n} = -\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^+}, \\
\boldsymbol{\sigma} \cdot \vec{n} = +\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^-},
\end{gather}
%
where $\vec{\lambda}$ is the vector of Lagrange multipliers corresponding to the tractions applied to the fault surface to generate the prescribed slip.

```{table} Mathematical notation for elasticity equation with infinitesimal strain and prescribed slip on faults.
:name: tab:notation:elasticity:prescribed:slip
| **Category**                   |   **Symbol**    | **Description**                                                                                   |
|:-------------------------------|:---------------:|:--------------------------------------------------------------------------------------------------|
| Unknowns                       |    $\vec{u}$    | Displacement field                                                                                |
|                                |    $\vec{v}$    | Velocity field                                                                                    |
|                                | $\vec{\lambda}$ | Lagrange multiplier field                                                                         |
| Derived quantities             |  $\boldsymbol{\sigma}$  | Cauchy stress tensor                                                                              |
|                                | $\boldsymbol{\epsilon}$ | Cauchy strain tensor                                                                              |
| Common constitutive parameters |     $\rho$      | Density                                                                                           |
|                                |      $\mu$      | Shear modulus                                                                                     |
|                                |       $K$       | Bulk modulus                                                                                      |
| Source terms                   |    $\vec{f}$    | Body force per unit volume, for example $\rho \vec{g}$                                            |
|                                |    $\vec{d}$    | Slip vector field on the fault corresponding to a jump in the displacement field across the fault |
```

:::{toctree}
quasistatic.md
dynamic.md
:::
