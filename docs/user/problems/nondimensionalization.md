# Nondimensionalization

Nondimensionalization permits solving a system equations across a vast range of scales.
We are primarily interested in elasticity-related boundary value problems, so we define displacement, length, rigidity, time, and temperature scales.

:::{seealso}
In the governing equations section, we discuss the nondimensionalization of each governing equation.

- [Elasticity](../governingeqns/elasticity/nondimensionalization.md)
- [Incompressible elasticity](../governingeqns/incompressible-elasticity/nondimensionalization.md)
- [Poroelasticity](../governingeqns/poroelasticity/nondimensionalization.md)
:::

## Scales

The Python `General` object is the generic interfaces to the underlying C++ `Scales` object that manages the scales used to nondimensionalize values.

:::{tip}
The default nondimensionalization is reasonable for many problems; however, it may be necessary to change the default values in some cases.
:::

## QuasistaticElasticity

`QuasistaticElasticity` provides convenient scales for nondimensionalizing quasistatic elasticity and incompressible elsticity boundary values problems.
The scales are specified in terms of the displacement scale, length scale, shear modulus, and time scale.

`displacement_scale`
: Nominal size of the displacment

`length_scale`
: Nominal size of the geometry controlling the deformation (for example, fault dimension or size of the domain)

`shear_modulus`
: Nominal shear modulus

`time_scale`
: Time scale of deformation (for example, viscoelastic relaxation time)

```{table} Nondimensional scales for quasistatic elasticity boundary value problems.
:name: tab:nondimensional:quasistatic:elasticity:scales
|  Scale   |      Name          | Default Value |
| :------: | :-------------:    | :-----------: |
|  length  | displacement_scale |    1 m        |
|  length  | length_scale       |  100 km       |
| pressure | shear_modulus      |   10 GPa      |
|   time   | time_scale         |  100 years    |
```

:::{admonition} Pyre User Interface
See [QuasistaticElasticity component](components/scales/QuasistaticElasticity.md).
:::

## DynamicElasticity

`DynamicElasticity` provides convenient scales for nondimensionalizing dynamic elasticity boundary values problems.
The scales are specified in terms of the displacement scale, length scale, shear modulus, and shear wave speed.
The time scale is computed based on the time it takes a shear wave ($v_s$)to propagate the length scale ($l_o$).

```{math}
t_o = l_o / v_s
```

`displacement_scale`
: Nominal size of the displacment

`length_scale`
: Nominal size of the geometry controlling the deformation (for example, fault dimension or size of the domain)

`shear_modulus`
: Nominal shear modulus

`shear_wave_speed`
: Nominal shear wave speed

```{table} Nondimensional scales for dynamic elasticity boundary value problems.
:name: tab:nondimensional:dynamic:elasticity:scales
|  Scale   |      Name          | Default Value |
| :------: | :-------------:    | :-----------: |
|  length  | displacement_scale |    1 m        |
|  length  | length_scale       |  100 km       |
| pressure | shear_modulus      |   10 GPa      |
|  --      | shear_wave_speed   |    3 km/s     |
```

:::{admonition} Pyre User Interface
See [DynamicElasticity component](components/scales/DynamicElasticity.md).
:::

## QuasistaticPoroelasticity

`QuasistaticPoroelasticity` provides convenient scales for nondimensionalizing quasistatic poroelasticity boundary values problems.
The scales are specified in terms of the displacement scale, length scale, shear modulus, viscosity, and permeability.
The time scale is computed based on the length scale ($l_o$), rigidity ($\mu_o$), viscosity ($\mu_{f_o}$), and permability ($k_o$).

```{math}
t_o = \frac{\mu_{f_o}}{k_o} \frac{l_o^2}{\mu_o}
```

`displacement_scale`
: Nominal size of the displacment

`length_scale`
: Nominal size of the geometry controlling the deformation (for example, fault dimension or size of the domain)

`shear_modulus`
: Nominal shear modulus

`viscosity`
: Nominal viscosity

`permeability`
: Nominal permeability

```{table} Nondimensional scales for quasistatic poroelasticity boundary value problems.
:name: tab:nondimensional:quasistatic:poroelasticity:scales
|  Scale   |      Name          | Default Value  |
| :------: | :-------------:    | :-----------:  |
|  length  | displacement_scale |    1 m         |
|  length  | length_scale       |  100 km        |
| pressure | shear_modulus      |   10 GPa       |
|  --      | viscosity          | 0.001 Pa-s     |
|  --      | permeability       | 1.0e-13 m$^2$  |
```

:::{admonition} Pyre User Interface
See [QuasistaticPoroelasticity component](components/scales/QuasistaticPoroelasticity.md).
:::
