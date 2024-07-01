# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`.
The `pylithapp.cfg` file contains numerous comments, so we only summarize the parameters here.

## Metadata, Mesh, and Output

The `pylithapp.metadata` section specifies metadata common to all simulations in the directory.
We control the verbosity of the output written to stdout using `journal.info`.
We set the parameters for importing the finite-element mesh from Gmsh in `pylithapp.mesh_generator`. 

## Physics

These quasi-static simulations solve the poroelasticity equation, so we have a solution field with  displacement, pressure, and volumetric strain subfields.
%
\begin{gather}
\vec{s} = \left(\vec{u} \quad p \quad \epsilon_v\right)^T, \\
\nabla \cdot \boldsymbol{\sigma}(\vec{u},p) = \vec{0}, \\
\frac{\partial \zeta(\vec{u},p)}{\partial t} + \nabla \cdot \vec{q}(p) = 0, \\
\nabla \cdot \vec{u} - \epsilon_{v} = 0.
\end{gather}

We specify a basis order of 2 for the displacement subfield, and 1 for each of the remaining subfields.
We also adjust the scales for nondimensionalization.

```{code-block} cfg
---
caption: Solution and nondimensionalization parameters for magma-2d examples.
---
[pylithapp.problem]
solution = pylith.problems.SolnDispPresTracStrain
defaults.quadrature_order = 2

[pylithapp.problem.solution.subfields]
displacement.basis_order = 2
pressure.basis_order = 1
trace_strain.basis_order = 1

[pylithapp.problem]
normalizer = spatialdata.units.NondimElasticQuasistatic
normalizer.length_scale = 100.0*m
normalizer.relaxation_time = 0.2*year
normalizer.shear_modulus = 10.0*GPa
```

We use the material properties in all of the simulations in this directory, so we specify them in `pylithapp.cfg` to avoid repeating the information in the file with parameters for each simulation.
We create an array of 2 materials; with uniform material properties we use `UniformDB` spatial databases and set the basis order of the auxiliary subfields to 0.

```{code-block} cfg
---
caption: Material parameters for the magma-2d examples.
---
[pylithapp.problem]
# We have two different poroelastic materials each with a linear bulk rheology.
materials = [crust, intrusion]
materials.crust = pylith.materials.Poroelasticity
materials.intrusion = pylith.materials.Poroelasticity

[pylithapp.problem.materials]
crust.bulk_rheology = pylith.materials.IsotropicLinearPoroelasticity
intrusion.bulk_rheology = pylith.materials.IsotropicLinearPoroelasticity

[pylithapp.problem.materials.crust]

# `label_value` must match the blocks in `bc.jou` Cubit Journal file.
description = crust
label_value = 1

# We will use uniform material properties, so we use the UniformDB
# spatial database.
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Poroelastic properties for the crust
db_auxiliary_field.values = [solid_density, fluid_density, fluid_viscosity, porosity, shear_modulus, drained_bulk_modulus, biot_coefficient, fluid_bulk_modulus, solid_bulk_modulus, isotropic_permeability]
db_auxiliary_field.data = [ 2500*kg/m**3, 1000*kg/m**3, 0.001*Pa*s, 0.01, 6.0*GPa, 10.0*GPa, 1.0, 2.0*GPa, 20.0*GPa, 1e-15*m**2]

# Set basis order to 0 for uniform properties and a basis order of 1 for Cauchy stress and strain.
auxiliary_subfields.body_force.basis_order = 0
auxiliary_subfields.solid_density.basis_order = 0
auxiliary_subfields.fluid_density.basis_order = 0
auxiliary_subfields.fluid_viscosity.basis_order = 0
auxiliary_subfields.gravitational_acceleration.basis_order = 0
auxiliary_subfields.porosity.basis_order = 0
derived_subfields.cauchy_strain.basis_order = 1
derived_subfields.cauchy_stress.basis_order = 1

[pylithapp.problem.materials.crust.bulk_rheology]
auxiliary_subfields.drained_bulk_modulus.basis_order = 0
auxiliary_subfields.shear_modulus.basis_order = 0
auxiliary_subfields.biot_coefficient.basis_order = 0
auxiliary_subfields.biot_modulus.basis_order = 0
auxiliary_subfields.isotropic_permeability.basis_order = 0


[pylithapp.problem.materials.intrusion]
# `label_value` must match the blocks in `bc.jou` Cubit Journal file.
description = Intrusion
label_value = 2

# We will use uniform material properties, so we use the UniformDB
# spatial database.
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Poroelastic properties
db_auxiliary_field.values = [solid_density, fluid_density, fluid_viscosity, porosity, shear_modulus, drained_bulk_modulus, biot_coefficient, fluid_bulk_modulus, solid_bulk_modulus, isotropic_permeability]
db_auxiliary_field.data = [ 2500*kg/m**3, 1000*kg/m**3, 0.001*Pa*s, 0.1, 6.0*GPa, 10.0*GPa, 0.8,  2.0*GPa, 20.0*GPa, 1e-13*m**2]

auxiliary_subfields.body_force.basis_order = 0
auxiliary_subfields.solid_density.basis_order = 0
auxiliary_subfields.fluid_density.basis_order = 0
auxiliary_subfields.fluid_viscosity.basis_order = 0
auxiliary_subfields.gravitational_acceleration.basis_order = 0
auxiliary_subfields.porosity.basis_order = 0
derived_subfields.cauchy_strain.basis_order = 1
derived_subfields.cauchy_stress.basis_order = 1

[pylithapp.problem.materials.intrusion.bulk_rheology]
# Set basis order to 0 for uniform properties
auxiliary_subfields.drained_bulk_modulus.basis_order = 0
auxiliary_subfields.shear_modulus.basis_order = 0
auxiliary_subfields.biot_coefficient.basis_order = 0
auxiliary_subfields.biot_modulus.basis_order = 0
auxiliary_subfields.isotropic_permeability.basis_order = 0
```
