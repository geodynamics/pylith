(sec:materials)=
# Materials (**`materials`**)

The material objects encapsulate the bulk behavior of the domain.
This includes both the governing equations and the associated bulk rheology.

## Specifying Material Properties

Associating material properties with a given cell involves several steps.

1.  In the mesh generation process, assign a material identifier to each cell.

2.  Define material property groups corresponding to each material identifier.
    In CUBIT/Trelis this is done by creating the blocks as part of the boundary conditions.

3.  Provide the settings for each material group in the parameters, i.e., `cfg` file.

4.  Specify the parameters for the material properties, e.g., linear variation in density with depth, using a spatial database file.
This allows variation of the material properties across cells with the same material identifier.

## Elasticity (`Elasticity`)

The *Elasticity* is used to solve the elasticity equation with or without inertia.
Whether inertia or body forces are included is determined by the *Elasticity* property settings.
Gravitational body forces are included if the **gravity_field** is set in the *Problem*.

:::{seealso}
[`Elasticity` Component](../components/materials/Elasticity.md)
:::

{numref}`tab:elasticity:rheologies` lists the bulk rheologies implemented for the elaticity equation.

```{table} Elasticity bulk rheologies.
:name: tab:elasticity:rheologies
|        Bulk Rheology        |  Description                                   |
|:----|:-----------------------------------------------|
| *IsotropicLinearElasticity* | Isotropic, linear elasticity                   |
|  *IsotropicLinearMaxwell*   | Isotropic, linear Maxwell viscoelasticity      |
| *IsotropicLinearGenMaxwell* | Isotropic, generalized Maxwell viscoelasticity |
|     *IsotropicPowerLaw*     | Isotropic, power-law viscoelasticity           |
|  *IsotropicDruckerPrager*   | Isotropic, Drucker-Prager elastoplasticity     |
```

:::{warning}
The *IsotropicDruckerPrager* rheology has not been implemented in this beta release.
:::

```{table} Auxiliary subfields for elasticity bulk rheologies.
:name: tab:elasticity:auxiliary:subfields
|       Subfield             |  L  | LM  | GM  |  PL | DP  |   Components           |
|:---------------------------|:---:|:---:|:---:|:---:|:---:|:-----------------------|
|                            |     |     |     |     |     |                        |
|                            |     |     |     |     |     |                        |
| density                    |  X  |  X  |  X  |  X  |  X  |                        |
| vp (P-wave speed)          |  X  |  X  |  X  |  X  |  X  |                        |
| vs (S-wave speed)          |  X  |  X  |  X  |  X  |  X  |                        |
| body_force                 |  O  |  O  |  O  |  O  |  O  | x, y, z                |
| gravitational_acceleration |  O  |  O  |  O  |  O  |  O  | x, y, z                |
| shear_modulus              |  I  |  I  |  I  |  I  |  I  |                        |
| bulk_modulus               |  I  |  I  |  I  |  I  |  I  |                        |
| reference_stress           |  O  |  O  |  O  |  O  |  O  | xx, yy, zz, xy, yz, xz |
| reference_strain           |  O  |  O  |  O  |  O  |  O  | xx, yy, zz, xy, yz, xz |
| maxwell_time               |     |  I  |     |     |     |                        |
| viscosity                  |     |  X  |     |     |     |                        |
| viscosity_1                |     |     |  X  |     |     |                        |
| viscosity_2                |     |     |  X  |     |     |                        |
| viscosity_3                |     |     |  X  |     |     |                        |
| shear_ratio_1              |     |     |  X  |     |     |                        |
| shear_ratio_2              |     |     |  X  |     |     |                        |
| shear_ratio_3              |     |     |  X  |     |     |                        |
| total_strain               |     |  X  |     |  X  |     | xx, yy, zz, xy, yz, xz |
| viscous_strain             |     |  X  |     |  X  |     | xx, yy, zz, xy, yz, xz |
| viscous_strain_1           |     |     |  X  |     |     | xx, yy, zz, xy, yz, xz |
| viscous_strain_2           |     |     |  X  |     |     | xx, yy, zz, xy, yz, xz |
| viscous_strain_3           |     |     |  X  |     |     | xx, yy, zz, xy, yz, xz |
| power_law_exponent         |     |     |     |  X  |     |                        |
| reference_strain_rate      |     |     |     |  X  |     |                        |
| reference_stress           |     |     |     |  X  |     |                        |
| cohesion                   |     |     |     |     |  X  |                        |
| friction_angle             |     |     |     |     |  X  |                        |
| dilatation_angle           |     |     |     |     |  X  |                        |
| alpha_yield                |     |     |     |     |  I  |                        |
| alpha_flow                 |     |     |     |     |  I  |                        |
| beta                       |     |     |     |     |  I  |                        |
| plastic_strain             |     |     |     |     |  X  | xx, yy, zz, xy, yz, xz |
```
X: required value in (db_auxiliary_fields) spatial database  
O: optional value in (db_auxiliary_fields) spatial database  
I: internal; computed from inputs  
L: isotropic, linear elasticity  
ML: isotropic linear Maxwell viscoelasticity  
GM: isotropic generalized linear Maxwell viscoelasticity  
PL: isotropic power-law viscoelasticity  
DP: isotropic Drucker-Prager elastoplasticity


```{table} Derived subfields for elasticity bulk rheologies.
:name: tab:elasticity:derived:subfields
|      Subfield | L   | LM  |  GM |  PL |  DP |  Components            |
|:--------------|:---:|:---:|:---:|:---:|:---:|:-----------------------|
| cauchy_stress |  ✓  |  ✓ |  ✓  |  ✓  |  ✓ | xx, yy, zz, xy, yz, xz |
| cauchy_strain |  ✓  |  ✓ |  ✓  |  ✓  |  ✓ | xx, yy, zz, yz, yz, xz |
```

## Incompressible Elasticity (`IncompressibleElasticity`)

Estimating realistic distributions of initial stress fields consistent with gravitational body forces can be quite difficult due to our lack of knowledge of the deformation history.
A simple way to approximate the lithostatic load is to solve for the stress field imposed by gravitational body forces assuming an incompressible elastic material.
This limits the volumetric deformation.
In this context we do not include inertia, so the *IncompressibleElasticity* object does not include an inertial term.
Gravitational body forces are included if the **gravity_field** is set in the *Problem*.

:::{seealso}
[`IncompressibleElasticity` Component](../components/materials/IncompressibleElasticity.md)
:::

{numref}`tab:incompressible:elasticity:rheologies` lists the bulk rheologies implemented for the elaticity equation.

```{table} Incompressible elasticity bulk rheologies.
:name: tab:incompressible:elasticity:rheologies
|   Bulk Rheology  | Description                                            |
|:----|:--------------------------------------------|
|  *IsotropicLinearincompElasticity*   | Isotropic, linear incompressible elasticity |
```

```{table} Auxiliary subfields for incompressible elasticity bulk rheologies.
:name: tab:incompressible:elasticity:auxiliary:subfields
|           Subfield         |   L | LM  |  GM |  PL |    Components          |
|:---------------------------|:---:|:---:|:---:|:---:|:-----------------------|
| density                    |  X  |  X  |  X  |  X  |                        |
| vp (P-wave speed)          |  X  |  X  |  X  |  X  |                        |
| vs (S-wave speed)          |  X  |  X  |  X  |  X  |                        |
| body_force                 |  O  |  O  |  O  |  O  | x, y, z                |
| gravitational_acceleration |  O  |  O  |  O  |  O  | x, y, z                |
| shear_modulus              |  I  |  I  |  I  |  I  |                        |
| bulk_modulus               |  I  |  I  |  I  |  I  |                        |
| reference_stress           |  O  |  O  |  O  |  O  | xx, yy, zz, xy, yz, xz |
| reference_strain           |  O  |  O  |  O  |  O  | xx, yy, zz, xy, yz, xz |
```
X: required value in (db_auxiliary_fields) spatial database  
O: optional value in (db_auxiliary_fields) spatial database  
I: internal; computed from inputs  
L: isotropic, linear elasticity  
ML: isotropic linear Maxwell viscoelasticity  
GM: isotropic generalized linear Maxwell viscoelasticity  
PL: isotropic power-law viscoelasticity  

```{table} Derived subfields for incompressible elasticity bulk rheologies.
:name: tab:incompressible:elasticity:derived:subfields
|      Subfield | L   | LM  |  GM |  PL | Components            |
|:--------------|:---:|:---:|:---:|:---:|:----------------------|
| cauchy_stress |  ✓  |  ✓ |  ✓  |  ✓  |xx, yy, zz, xy, yz, xz |
| cauchy_strain |  ✓  |  ✓ |  ✓  |  ✓  |xx, yy, zz, yz, yz, xz |
```
