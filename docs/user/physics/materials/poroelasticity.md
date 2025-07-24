(sec-user-physics-poroelasticity)=
# Poroelasticity

You can use the `Poroelasticity` component to solve the poroelasticity equation with or without inertia.
Whether inertia or body forces are included is determined by the `Poroelasticity` property settings.
Gravitational body forces are included if the `gravity_field` is set in the `Problem`.
{numref}`tab:poroelasticity:rheologies` lists the poroelastic bulk rheology implemented for the poroelaticity equation.

```{table} Elasticity bulk rheologies.
:name: tab:poroelasticity:rheologies
| Bulk Rheology                   | Description                      |
| :------------------------------ | :------------------------------- |
| `IsotropicLinearPoroelasticity` | Isotropic, linear poroelasticity |
```

```{table} Properties defining elasticity bulk rheologies.
:name: tab:poroelasticity:auxiliary:subfields
| Subfield                     |   L   | Components             |
| :--------------------------- | :---: | :--------------------- |
| `solid_density`              |   X   |                        |
| `fluid_density`              |   X   |                        |
| `fluid_viscosity`            |   X   |                        |
| `porosity`                   |   X   |                        |
| `body_force`                 |   O   | x, y, z                |
| `gravitational_acceleration` |   O   | x, y, z                |
| `shear_modulus`              |   X   |                        |
| `drained_bulk_modulus`       |   X   |                        |
| `fluid_bulk_modulus`         |   X   |                        |
| `biot_coefficient`           |   X   |                        |
| `isotropic_permeability`     |   X   |                        |
| `tensor_permeability`        |   O   | xx, yy, zz, xy, yz, xz |
| `reference_stress`           |   O   | xx, yy, zz, xy, yz, xz |
| `reference_strain`           |   O   | xx, yy, zz, xy, yz, xz |
```

X: required value in auxiliary field spatial database  
O: optional value in auxiliary field spatial database  
L: isotropic, linear poroelasticity  

```{table} Derived subfields that are available for output for poroelasticity bulk rheologies.
:name: tab:poroelasticity:derived:subfields
| Subfield        |   L   | Components             |
| :-------------- | :---: | :--------------------- |
| `cauchy_stress` |   X   | xx, yy, zz, xy, yz, xz |
| `cauchy_strain` |   X   | xx, yy, zz, xy, yz, xz |
| `bulk_density`  |   X   |                        |
| `water_content` |   X   |                        |
```

When porosity is enabled as a state variable, it will be included in the output along with the derived subfields.

:::{seealso}
See [`Poroelasticity` Component](../../components/materials/Poroelasticity.md) for the Pyre properties and facilities and configuration examples.
:::
