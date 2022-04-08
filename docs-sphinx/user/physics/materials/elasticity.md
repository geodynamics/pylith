(sec-user-physics-elasticity)=
# Elasticity

You can use the `Elasticity` component to solve the elasticity equation with or without inertia.
Whether inertia or body forces are included is determined by the `Elasticity` property settings.
Gravitational body forces are included if the `gravity_field` is set in the `Problem`.
{numref}`tab:elasticity:rheologies` lists the various elastic, viscoelastic, and elastoplastic bulk constitutive models implemented for the elaticity equation.

```{table} Elasticity bulk rheologies.
:name: tab:elasticity:rheologies
| Bulk Rheology               |  Description                                   |
|:----------------------------|:-----------------------------------------------|
| `IsotropicLinearElasticity` | Isotropic, linear elasticity                   |
| `IsotropicLinearMaxwell`    | Isotropic, linear Maxwell viscoelasticity      |
| `IsotropicLinearGenMaxwell` | Isotropic, generalized Maxwell viscoelasticity |
| `IsotropicPowerLaw`         | Isotropic, power-law viscoelasticity           |
| `IsotropicDruckerPrager`    | Isotropic, Drucker-Prager elastoplasticity     |
```

:::{warning}
The `IsotropicDruckerPrager` rheology has not yet been implemented in PyLith v3.
:::

```{table} Properties defining elasticity bulk rheologies.
:name: tab:elasticity:auxiliary:subfields
|       Subfield               |  L  | LM  | GM  |  PL | DP  |   Components           |
|:-----------------------------|:---:|:---:|:---:|:---:|:---:|:-----------------------|
| `density`                    |  X  |  X  |  X  |  X  |  X  |                        |
| `vp` (P-wave speed)          |  X  |  X  |  X  |  X  |  X  |                        |
| `vs` (S-wave speed)          |  X  |  X  |  X  |  X  |  X  |                        |
| `body_force`                 |  O  |  O  |  O  |  O  |  O  | x, y, z                |
| `gravitational_acceleration` |  O  |  O  |  O  |  O  |  O  | x, y, z                |
| `shear_modulus`              |  I  |  I  |  I  |  I  |  I  |                        |
| `bulk_modulus`               |  I  |  I  |  I  |  I  |  I  |                        |
| `reference_stress`           |  O  |  O  |  O  |  O  |  O  | xx, yy, zz, xy, yz, xz |
| `reference_strain`           |  O  |  O  |  O  |  O  |  O  | xx, yy, zz, xy, yz, xz |
| `maxwell_time`               |     |  I  |     |     |     |                        |
| `viscosity`                  |     |  X  |     |     |     |                        |
| `viscosity_1`                |     |     |  X  |     |     |                        |
| `viscosity_2`                |     |     |  X  |     |     |                        |
| `viscosity_3`                |     |     |  X  |     |     |                        |
| `shear_ratio_1`              |     |     |  X  |     |     |                        |
| `shear_ratio_2`              |     |     |  X  |     |     |                        |
| `shear_ratio_3`              |     |     |  X  |     |     |                        |
| `total_strain`               |     |  X  |     |  X  |     | xx, yy, zz, xy, yz, xz |
| `viscous_strain`             |     |  X  |     |  X  |     | xx, yy, zz, xy, yz, xz |
| `viscous_strain_1`           |     |     |  X  |     |     | xx, yy, zz, xy, yz, xz |
| `viscous_strain_2`           |     |     |  X  |     |     | xx, yy, zz, xy, yz, xz |
| `viscous_strain_3`           |     |     |  X  |     |     | xx, yy, zz, xy, yz, xz |
| `power_law_exponent`         |     |     |     |  X  |     |                        |
| `reference_strain_rate`      |     |     |     |  X  |     |                        |
| `reference_stress`           |     |     |     |  X  |     |                        |
| `cohesion`                   |     |     |     |     |  X  |                        |
| `friction_angle`             |     |     |     |     |  X  |                        |
| `dilatation_angle`           |     |     |     |     |  X  |                        |
| `alpha_yield`                |     |     |     |     |  I  |                        |
| `alpha_flow`                 |     |     |     |     |  I  |                        |
| `beta`                       |     |     |     |     |  I  |                        |
| `plastic_strain`             |     |     |     |     |  X  | xx, yy, zz, xy, yz, xz |
```

X: required value in auxiliary field spatial database  
O: optional value in auxiliary field spatial database  
I: internal auxiliary subfield; computed from spatial database values
L: isotropic, linear elasticity  
ML: isotropic linear Maxwell viscoelasticity  
GM: isotropic generalized linear Maxwell viscoelasticity  
PL: isotropic power-law viscoelasticity  
DP: isotropic Drucker-Prager elastoplasticity

```{table} Derived subfields that are available for output for elasticity bulk rheologies.
:name: tab:elasticity:derived:subfields
|      Subfield | L   | LM  |  GM |  PL |  DP |  Components            |
|:--------------|:---:|:---:|:---:|:---:|:---:|:-----------------------|
| `cauchy_stress` |  ✓  |  ✓ |  ✓  |  ✓  |  ✓ | xx, yy, zz, xy, yz, xz |
| `cauchy_strain` |  ✓  |  ✓ |  ✓  |  ✓  |  ✓ | xx, yy, zz, xy, yz, xz |
```

:::{seealso}
See [`Elasticity` Component](../../components/materials/Elasticity.md) for the Pyre properties and facilities and configuration examples.
:::
