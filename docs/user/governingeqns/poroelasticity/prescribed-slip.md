# Infinitesimal Strain and Prescribed Fault Slip

We apply the prescribed fault slip formulation that we used for elasticity to poroelasticity.
We add a boundary condition to the poroelasticity equation prescribing the jump in the displacement field across the fault,
%
```{math}
\begin{gathered}
  -\vec{u}^+ + \vec{u}^- + \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f,
\end{gathered}
```
%
where $\vec{u}^+$ is the displacement vector on the "positive" side of the fault, $\vec{u}^-$ is the displacement vector on the "negative" side of the fault, $\vec{d}$ is the slip vector on the fault, and $\vec{n}$ is the fault normal which points from the negative side of the fault to the positive side of the fault.
We enforce the jump in displacements across the fault using a Lagrange multiplier corresponding to equal and opposite tractions on the two sides of the fault.

:::{warning}
In this formulation, the fault acts as a barrier to fluid flow.
That is, there is no coupling in fluid pressure across the fault.
:::

```{table} Mathematical notation for poroelasticity with infinitesimal strain and prescribed fault slip.
:name: tab:notation:poroelasticity:prescribed:slip

| **Category**                   |       **Symbol**        | **Description**                                                                                               |
| :----------------------------- | :---------------------: | :------------------------------------------------------------------------------------------------------------ |
| Unknowns                       |        $\vec{u}$        | Displacement field                                                                                            |
|                                |        $\vec{v}$        | Velocity field                                                                                                |
|                                |           $p$           | Pressure field (corresponds to pore fluid pressure)                                                           |
|                                |     $\epsilon_{v}$      | Volumetric (trace) strain                                                                                     |
|                                |     $\vec{\lambda}$     | Lagrange multiplier for fault                                                                                 |
| Derived quantities             |  $\boldsymbol{\sigma}$  | Cauchy stress tensor                                                                                          |
|                                | $\boldsymbol{\epsilon}$ | Cauchy strain tensor                                                                                          |
|                                |         $\zeta$         | Variation of fluid content (variation of fluid vol. per unit vol. of PM), $\alpha \epsilon_{v} + \frac{p}{M}$ |
|                                |       $\rho_{b}$        | Bulk density, $\left(1 - \phi\right) \rho_{s} + \phi \rho_{f}$                                                |
|                                |        $\vec{q}$        | Darcy flux, $-\frac{\boldsymbol{k}}{\mu_{f}} \cdot \left(\nabla p - \vec{f}_{f}\right)$                       |
| Common constitutive parameters |       $\rho_{s}$        | Solid (matrix) density                                                                                        |
|                                |       $\rho_{f}$        | Fluid density                                                                                                 |
|                                |        $\mu_{f}$        | Fluid viscosity                                                                                               |
|                                |         $\phi$          | Porosity                                                                                                      |
|                                |          $\mu$          | Shear modulus                                                                                                 |
|                                |         $K_{d}$         | Drained bulk modulus                                                                                          |
|                                |        $\alpha$         | Biot coefficient, $1 - \frac{K_{d}}{K_{s}}$                                                                   |
|                                |           $M$           | Biot modulus                                                                                                  |
|                                |    $\boldsymbol{k}$     | Permeability                                                                                                  |
| Source terms                   |        $\vec{f}$        | Body force per unit volume, for example: $\rho_{b} \vec{g}$                                                   |
|                                |      $\vec{f}_{f}$      | Fluid body force, for example: $\rho_{f} \vec{g}$                                                             |
|                                |        $\gamma$         | Source density; rate of injected fluid per unit volume of the porous solid                                    |
```

:::{toctree}
prescribed-slip-quasistatic.md
:::
