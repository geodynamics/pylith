# QuasistaticElasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.scales.QuasistaticElasticity`
:Journal name: `quasistaticelasticity`

Convenience object for nondimensionalizing quasi-static elasticity problems.

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
* `time_scale`=\<dimensional\>: Time scale of boundary value problem (for example, viscoelastic relaxation time).
  - **default value**: 3.15576e+09*s
  - **current value**: 3.15576e+09*s, from {default}
  - **validator**: (greater than 0*s)

## Example

Example of setting `QuasistaticElasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[normalizer]
length_scale = 100.0*km
displacement_scale = 50.0*km
shear_modulus = 25.0*GPa
time_scale = 100.0*year
:::

