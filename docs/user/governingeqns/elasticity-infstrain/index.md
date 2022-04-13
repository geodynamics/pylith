# Elasticity with Infinitesimal Strain and No Faults

We begin with the elasticity equation including the inertial term,

```{math}
:label: eqn:elasticity:strong:form
\rho \frac{\partial^2\vec{u}}{\partial t^2} - \vec{f}(\vec{x},t) - \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} (\vec{u}) = \vec{0} \text{ in }\Omega,
```

```{math}
:label: eqn:bc:Neumann
\boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on }\Gamma_\tau,
```

```{math}
:label: eqn:bc:Dirichlet
\vec{u} = \vec{u}_0(\vec{x},t) \text{ on }\Gamma_u,
```

where $\vec{u}$ is the displacement vector, $\rho$ is the mass density, $\vec{f}$ is the body force vector, $\boldsymbol{\sigma}$ is the Cauchy stress tensor, $\vec{x}$ is the spatial coordinate, and $t$ is time. We specify tractions $\vec{\tau}$ on boundary $\Gamma_\tau$, and displacements $\vec{u}_0$ on boundary $\Gamma_u$.
Because both $\vec{\tau}$ and $\vec{u}$ are vector quantities, there can be some spatial overlap of boundaries $\Gamma_\tau$ and $\Gamma_u$; however, a degree of freedom at any location cannot be associated with both prescribed displacements (Dirichlet) and traction (Neumann) boundary conditions simultaneously.

```{table} Mathematical notation for elasticity equation with infinitesimal strain.
:name: tab:notation:elasticity

| **Category**                   |   **Symbol**    | **Description**                                        |
|:-------------------------------|:---------------:|:-------------------------------------------------------|
| Unknowns                       |    $\vec{u}$    | Displacement field                                     |
|                                |    $\vec{v}$    | Velocity field                                         |
| Derived quantities             |  $\boldsymbol{\sigma}$  | Cauchy stress tensor                                   |
|                                | $\boldsymbol{\epsilon}$ | Cauchy strain tensor                                   |
| Common constitutive parameters |     $\rho$      | Density                                                |
|                                |      $\mu$      | Shear modulus                                          |
|                                |       $K$       | Bulk modulus                                           |
| Source terms                   |    $\vec{f}$    | Body force per unit volume, for example $\rho \vec{g}$ |
```

:::{toctree}
quasistatic.md
dynamic.md
:::
