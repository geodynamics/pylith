# QuasistaticPoroelasticity

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.scales.QuasistaticPoroelasticity`
:Journal name: `quasistaticporoelasticity`

Convenience object for nondimensionalizing quasi-static poroelasticity problems.

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
* `permeability`=\<dimensional\>: Nominal permeability in boundary value problem.
  - **default value**: 1e-13*m**2
  - **current value**: 1e-13*m**2, from {default}
  - **validator**: (greater than 0*m**-1*kg*s**-2)
* `shear_modulus`=\<dimensional\>: Nominal shear modulus in boundary value problem.
  - **default value**: 1e+10*m**-1*kg*s**-2
  - **current value**: 1e+10*m**-1*kg*s**-2, from {default}
  - **validator**: (greater than 0*m**-1*kg*s**-2)
* `viscosity`=\<dimensional\>: Nominal fluid viscosity in boundary value problem.
  - **default value**: 0.001*m**-1*kg*s**-1
  - **current value**: 0.001*m**-1*kg*s**-1, from {default}
  - **validator**: (greater than 0*m**-1*kg*s**-2)

## Example

Example of setting `QuasistaticPoroelasticity` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[normalizer]
length_scale = 100.0*km
displacement_scale = 50.0*km
shear_modulus = 25.0*GPa
viscosity = 0.001*Pa*s
permeability = 1.0e-12*m**2
:::

