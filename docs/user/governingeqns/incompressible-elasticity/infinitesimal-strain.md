# Infinitesimal Strain (Bathe) and No Faults

In this section we apply a similar approach to the one we use for the elasticity equation to the case of an incompressible material.
As the bulk modulus ($K$) approaches infinity, the volumetric strain ($\mathop{\mathrm{Tr}}(\epsilon)$) approaches zero and the pressure remains finite, $p = -K \mathop{\mathrm{Tr}}(\epsilon)$.
We consider pressure $p$ as an independent variable and decompose the stress into the pressure and deviatoric components.
As a result, we write the stress tensor in terms of both the displacement and pressure fields,
%
\begin{equation}
\boldsymbol{\sigma}(\vec{u},p) = \boldsymbol{\sigma}^\mathit{dev}(\vec{u}) - p\boldsymbol{I}.
\end{equation}

The strong form is
%
\begin{gather}
  % Solution
  \vec{s}^T = \left( \vec{u} \quad \ p \right)^T, \\
  % Elasticity
\rho \frac{\partial^2\vec{u}}{\partial t^2} - \vec{f}(\vec{x},t) - \left(\boldsymbol{\sigma}^\mathit{dev}(\vec{u}) - p\boldsymbol{I}\right) = \vec{0} \text{ in }\Omega, \\
  % Pressure
  \vec{\nabla} \cdot \vec{u} + \frac{p}{K} = 0 \text{ in }\Omega, \\
  % Neumann
  \boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau} \text{ on }\Gamma_\tau, \\
  % Dirichlet
  \vec{u} = \vec{u}_0 \text{ on }\Gamma_u, \\
  p = p_0 \text{ on }\Gamma_p.
\end{gather}
%

```{table} Mathematical notation for incompressible elasticity with infinitesimal strain
:name: tab:notation:incompressible:elasticity

| Category                       |               Symbol               | Description                                                |
| :----------------------------- | :--------------------------------: | :--------------------------------------------------------- |
| Unknowns                       |             $\vec{u}$              | Displacement field                                         |
|                                |                $p$                 | Pressure field ($p>0$ corresponds to negative mean stress) |
| Derived quantities             |       $\boldsymbol{\sigma}$        | Cauchy stress tensor                                       |
|                                | $\boldsymbol{\sigma}^\mathit{dev}$ | Cauchy deviatoric stress tensor                            |
|                                |      $\boldsymbol{\epsilon}$       | Cauchy strain tensor                                       |
| Common constitutive parameters |               $\rho$               | Density                                                    |
|                                |               $\mu$                | Shear modulus                                              |
|                                |                $K$                 | Bulk modulus                                               |
| Source terms                   |             $\vec{f}$              | Body force per unit volume, for example $\rho \vec{g}$     |
```

:::{toctree}
infinitesimal-strain-quasistatic.md
:::
