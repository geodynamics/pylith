# Poroelasticity with Infinitesimal Strain and No Faults

We base this formulation for poroelasticity on Zheng et al. and Detournay and Cheng (1993).
We assume a slightly compressible fluid that completely saturates a porous solid, undergoing infinitesimal strain.

We begin with the conservation of linear momentum, including inertia, borrowed from linear elasticity:
%
\begin{equation}
  \rho_s\frac{\partial^2 \vec{u}}{\partial t^2} = \vec{f}(t) + \nabla \cdot \boldsymbol{\sigma}(\vec{u},p).
\end{equation}
%
Enforcing mass balance of the fluid gives
%
\begin{gather}
  \frac{\partial \zeta(\vec{u},p)}{\partial t} + \nabla \cdot \vec{q}(p) =
  \gamma(\vec{x},t) \text{ in } \Omega, \\
%
  \vec{q} \cdot \vec{n} = q_0(\vec{x},t) \text{ on }\Gamma_q, \\
%
  p = p_0(\vec{x},t) \text{ on }\Gamma_p,
\end{gather}
%
where $\zeta$ is the variation in fluid content, $\vec{q}$ is the rate of fluid volume crossing a unit area of the porous solid, $\gamma$ is the rate of injected fluid per unit volume of the porous solid, $q_0$ is the outward fluid velocity normal to the boundary $\Gamma_q$, and $p_0$ is the fluid pressure on boundary $\Gamma_p$.

We require the fluid flow to follow Darcy's law (Navier-Stokes equation neglecting inertial effects),
%
\begin{equation}
  \vec{q}(p) = -\frac{\boldsymbol{k}}{\mu_{f}}(\nabla p - \vec{f}_f),
\end{equation}
%
where $\boldsymbol{k}$ is the intrinsic permeability, $\mu_f$ is the viscosity of the fluid, $p$ is the fluid pressure, and $\vec{f}_f$ is the body force in the fluid.
If gravity is included in a problem, then usually $\vec{f}_f = \rho_f \vec{g}$, where $\rho_f$ is the density of the fluid and $\vec{g}$ is the gravitational acceleration vector.

## Constitutive Behavior

We assume linear elasticity for the solid phase, so the constitutive behavior can be expressed as
%
\begin{equation}
  \boldsymbol{\sigma}(\vec{u},p) = \boldsymbol{C} : \boldsymbol{\epsilon} - \alpha p \boldsymbol{I},
\end{equation}
%
where $\boldsymbol{\sigma}$ is the stress tensor, $\boldsymbol{C}$ is the tensor of elasticity constants, $\alpha$ is the Biot coefficient (effective stress coefficient), $\boldsymbol{\epsilon}$ is the strain tensor, and $\boldsymbol{I}$ is the identity tensor.
For this case, we will assume that the material properties are isotropic, resulting in the following formulation for the stress tensor:
%
\begin{equation}
  \boldsymbol{\sigma}(\vec{u},p) = \boldsymbol{C}:\boldsymbol{\epsilon} - \alpha p \boldsymbol{I}
  = \lambda \boldsymbol{I} \epsilon_{v} + 2 \mu - \alpha \boldsymbol{I} p
\end{equation}
%
where $\lambda$ and $\mu$ are Lam&eacute;'s parameters, $\lambda = K_{d} - \frac{2 \mu}{3}$, $\mu$ is the shear modulus, and the volumetric strain is defined as $\epsilon_{v} = \nabla \cdot \vec{u}$.

For the constitutive behavior of the fluid, we use the volumetric strain to couple the fluid-solid behavior,
%
\begin{gather}
  \zeta(\vec{u},p) = \alpha \mathop{\mathrm{Tr}}({\boldsymbol{\epsilon}}) + \frac{p}{M}, \\
%
  \frac{1}{M} = \frac{\alpha-\phi}{K_s} + \frac{\phi}{K_f},
\end{gather}
%
where $1/M$ is the specific storage coefficient at constant strain, $K_s$ is the bulk modulus of the solid, and $K_f$ is the bulk modulus of the fluid.
We can write the trace of the strain tensor as the dot product of the gradient and displacement field, so we have
%
\begin{equation}
  \zeta(\vec{u},p) = \alpha (\nabla \cdot \vec{u}) + \frac{p}{M}.
\end{equation}
%
```{table} Mathematical notation for poroelasticity with infinitesimal strain.
:name: tab:notation:poroelasticity

| **Category**                   |   **Symbol**    | **Description**                                                                                               |
|:-------------------------------|:---------------:|:--------------------------------------------------------------------------------------------------------------|
| Unknowns                       |    $\vec{u}$    | Displacement field                                                                                            |
|                                |    $\vec{v}$    | Velocity field                                                                                                |
|                                |       $p$       | Pressure field (corresponds to pore fluid pressure)                                                           |
|                                | $\epsilon_{v}$  | Volumetric (trace) strain                                                                                     |
| Derived quantities             |  $\boldsymbol{\sigma}$  | Cauchy stress tensor                                                                                  |
|                                | $\boldsymbol{\epsilon}$ | Cauchy strain tensor                                                                                  |
|                                |     $\zeta$     | Variation of fluid content (variation of fluid vol. per unit vol. of PM), $\alpha \epsilon_{v} + \frac{p}{M}$ |
|                                |   $\rho_{b}$    | Bulk density, $\left(1 - \phi\right) \rho_{s} + \phi \rho_{f}$                                                |
|                                |    $\vec{q}$    | Darcy flux, $-\frac{\boldsymbol{k}}{\mu_{f}} \cdot \left(\nabla p - \vec{f}_{f}\right)$                       |
|                                |       $M$       | Biot modulus                                                                                                  |
| Common constitutive parameters |   $\rho_{f}$    | Fluid density                                                                                                 |
|                                |   $\rho_{s}$    | Solid (matrix) density                                                                                        |
|                                |     $\phi$      | Porosity                                                                                                      |
|                                |    $\boldsymbol{k}$     | Permeability                                                                                          |
|                                |    $\mu_{f}$    | Fluid viscosity                                                                                               |
|                                |     $K_{s}$     | Solid grain bulk modulus                                                                                      |
|                                |     $K_{f}$     | Fluid bulk modulus                                                                                            |
|                                |     $K_{d}$     | Drained bulk modulus                                                                                          |
|                                |    $\alpha$     | Biot coefficient, $1 - \frac{K_{d}}{K_{s}}$                                                                   |
| Source terms                   |    $\vec{f}$    | Body force per unit volume, for example: $\rho_{b} \vec{g}$                                                   |
|                                |  $\vec{f}_{f}$  | Fluid body force, for example: $\rho_{f} \vec{g}$                                                             |
|                                |    $\gamma$     | Source density; rate of injected fluid per unit volume of the porous solid                                    |
```

:::{toctree}
quasistatic.md
dynamic.md
:::
