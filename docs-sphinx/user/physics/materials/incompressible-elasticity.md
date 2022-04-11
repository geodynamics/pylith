(sec-user-physics-incompressible-elasticity)=
# Incompressible Elasticity

You can use the `IncompressibleElasticity` component to solve the quasistatic incompressible elasticity equation.
Estimating realistic distributions of initial stress fields consistent with gravitational body forces can be quite difficult due to our lack of knowledge of the deformation history.
A simple way to approximate the lithostatic load is to solve for the stress field imposed by gravitational body forces assuming an incompressible elastic material.
This limits the volumetric deformation.
In this context we do not include inertia, so the `IncompressibleElasticity` component does not include an inertial term.
Gravitational body forces are included if the `gravity_field` is set in the `Problem`.
{numref}`tab:incompressible:elasticity:rheologies` lists the elastic bulk rheology implemented for the incompressible elaticity equation.

```{table} Incompressible elasticity bulk rheology.
:name: tab:incompressible:elasticity:rheologies
|   Bulk Rheology                   | Description                                 |
|:----------------------------------|:--------------------------------------------|
| `IsotropicLinearincompElasticity` | Isotropic, linear incompressible elasticity |
```

```{table} Properties defining incompressible elasticity bulk rheologies.
:name: tab:incompressible:elasticity:auxiliary:subfields
|           Subfield         |   L | LM  |  GM |  PL |    Components          |
|:---------------------------|:---:|:---:|:---:|:---:|:-----------------------|
| `density`                    |  X  |  X  |  X  |  X  |                        |
| `vp` (P-wave speed)          |  X  |  X  |  X  |  X  |                        |
| `vs` (S-wave speed)          |  X  |  X  |  X  |  X  |                        |
| `body_force`                 |  O  |  O  |  O  |  O  | x, y, z                |
| `gravitational_acceleration` |  O  |  O  |  O  |  O  | x, y, z                |
| `shear_modulus`              |  I  |  I  |  I  |  I  |                        |
| `bulk_modulus`               |  I  |  I  |  I  |  I  |                        |
| `reference_stress`           |  O  |  O  |  O  |  O  | xx, yy, zz, xy, yz, xz |
| `reference_strain`           |  O  |  O  |  O  |  O  | xx, yy, zz, xy, yz, xz |
```
X: required value in auxiliary field spatial database  
O: optional value in auxiliary field spatial database  
I: internal auxiliary subfield; computed from spatial database values
L: isotropic, linear elasticity  
ML: isotropic linear Maxwell viscoelasticity  
GM: isotropic generalized linear Maxwell viscoelasticity  
PL: isotropic power-law viscoelasticity  

```{table} Derived subfields that are available for output for incompressible elasticity bulk rheologies.
:name: tab:incompressible:elasticity:derived:subfields
|      Subfield | L   | LM  |  GM |  PL | Components            |
|:--------------|:---:|:---:|:---:|:---:|:----------------------|
| `cauchy_stress` |  ✓  |  ✓ |  ✓  |  ✓  |xx, yy, zz, xy, yz, xz |
| `cauchy_strain` |  ✓  |  ✓ |  ✓  |  ✓  |xx, yy, zz, xy, yz, xz |
```

:::{seealso}
See [`IncompressibleElasticity` Component](../../components/materials/IncompressibleElasticity.md) for the Pyre properties and facilities and configuration examples.
:::
