# DynamicElasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.scales.DynamicElasticity`
:Journal name: `dynamicelasticity`

Convenience object for nondimensionalizing dynamic elasticity problems.

Implements `General`.

## Pyre Properties

* `displacement_scale`=\<dimensional\>: Nominal displacement scale in boundary value problem.
  - **default value**: 1*m
  - **current value**: 1*m, from {default}
  - **validator**: (greater than 0*m)
* `length_scale`=\<dimensional\>: Length scale in boundary value problem (size of feature controlling displacement, fault).
  - **default value**: 100000*m
  - **current value**: 100000*m, from {default}
  - **validator**: (greater than 0*m)
* `shear_modulus`=\<dimensional\>: Nominal shear modulus in boundary value problem.
  - **default value**: 1e+10*m**-1*kg*s**-2
  - **current value**: 1e+10*m**-1*kg*s**-2, from {default}
  - **validator**: (greater than 0*m**-1*kg*s**-2)
* `shear_wave_speed`=\<dimensional\>: Nominal shear wave speed in boundary value problem.
  - **default value**: 3000*m*s**-1
  - **current value**: 3000*m*s**-1, from {default}
  - **validator**: (greater than 0*s)

## Example

Example of setting `DynamicElasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[normalizer]
length_scale = 1.0*km
displacement_scale = 50.0*km
shear_modulus = 25.0*GPa
velocity_scale = 3.0*km/s
:::

