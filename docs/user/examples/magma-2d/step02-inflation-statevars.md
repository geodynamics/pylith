# Step 2: Magma inflation with evolution of porosity

*New in v4.0.0*

```{include} step02_inflation_statevars-synopsis.md
```

## Simulation parameters

We extend the simulation in Step 1 by including evolution of the porosity, which depends on the time derivative of the pressure and trace strain.
We also compute the deformation relative to a uniform reference compressive pressure of 5 MPa to illustrate how to use a reference state with poroelasticity.
We use the same initial conditions and boundary conditions as in Step 1.

Because the evolution of porosity depends on the time derivative of the solution subfields, we need to include the time derivatives in the solution field.
As a result, we have 6 subfields in our solution field.

```{code-block} cfg
---
caption: Solution subfields for Step 2.
---
# Poroelasticity with porosity state variable requires solution with time derivatives 
solution = pylith.problems.SolnDispPresTracStrainVelPdotTdot

# Set basis order for all solution subfields
[pylithapp.problem.solution.subfields]
displacement.basis_order = 2
pressure.basis_order = 1
trace_strain.basis_order = 1
velocity.basis_order = 2
pressure_t.basis_order = 1
trace_strain_t.basis_order = 1
```

```{code-block} cfg
---
caption: Material parameters for poroelasticity with state variables and reference state for Step 2.
---
[pylithapp.problem.materials.crust]
use_state_variables = True

db_auxiliary_field.values = [
    solid_density, fluid_density, fluid_viscosity, porosity, shear_modulus, drained_bulk_modulus, biot_coefficient, fluid_bulk_modulus, solid_bulk_modulus, isotropic_permeability,
    reference_stress_xx, reference_stress_yy, reference_stress_zz, reference_stress_xy,
    reference_strain_xx, reference_strain_yy, reference_strain_zz, reference_strain_xy
    ]
db_auxiliary_field.data   = [
    2500*kg/m**3, 1000*kg/m**3, 0.001*Pa*s, 0.01, 6.0*GPa, 10.0*GPa, 1.0, 2.0*GPa, 20.0*GPa, 1e-15*m**2,
    -5.0*MPa, -5.0*MPa, -5.0*MPa, 0.0*MPa,
    0.0, 0.0, 0.0, 0.0
    ]

auxiliary_subfields.porosity.basis_order = 1

[pylithapp.problem.materials.crust.bulk_rheology]
use_reference_state = True


[pylithapp.problem.materials.intrusion]
use_state_variables = True

db_auxiliary_field.values = [
    solid_density, fluid_density, fluid_viscosity, porosity, shear_modulus, drained_bulk_modulus, biot_coefficient, fluid_bulk_modulus, solid_bulk_modulus, isotropic_permeability,
    reference_stress_xx, reference_stress_yy, reference_stress_zz, reference_stress_xy,
    reference_strain_xx, reference_strain_yy, reference_strain_zz, reference_strain_xy
    ]
db_auxiliary_field.data   = [
    2500*kg/m**3,  1000*kg/m**3, 0.001*Pa*s, 0.1, 6.0*GPa, 10.0*GPa, 0.8, 2.0*GPa, 20.0*GPa, 1e-13*m**2,
    -5.0*MPa, -5.0*MPa, -5.0*MPa, 0.0*Pa,
    0.0, 0.0, 0.0, 0.0
    ]

auxiliary_subfields.porosity.basis_order = 1

[pylithapp.problem.materials.intrusion.bulk_rheology]
use_reference_state = True
```

## Running the simulation

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02_inflation.cfg

# The output should look something like the following.

 >> /Users/baagaard/software/unix/py310-venv/pylith-debug/lib/python3.10/site-packages/pylith/apps/PyLithApp.py:84:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> /Users/baagaard/software/unix/py310-venv/pylith-debug/lib/python3.10/site-packages/pylith/meshio/MeshIOObj.py:43:read
 -- meshiocubit(info)
 -- Reading finite-element mesh
 >> /Users/baagaard/src/cig/pylith/libsrc/pylith/meshio/MeshIOCubit.cc:156:void pylith::meshio::MeshIOCubit::_readVertices(ExodusII &, scalar_array *, int *, int *) const
 -- meshiocubit(info)
 -- Component 'reader': Reading 747 vertices.

# -- many lines omitted --

 >> /Users/baagaard/software/unix/py310-venv/pylith-debug/lib/python3.10/site-packages/pylith/problems/TimeDependent.py:137:run
 -- timedependent(info)
 -- Solving problem.
0 TS dt 1. time -1.
    0 SNES Function norm 7.521665654021e-01
    Linear solve converged due to CONVERGED_ATOL iterations 1
    1 SNES Function norm 5.516056137722e-02
    Linear solve converged due to CONVERGED_ATOL iterations 1
    2 SNES Function norm 1.166275103468e-02
    Linear solve converged due to CONVERGED_ATOL iterations 1
    3 SNES Function norm 3.829709526109e-03
    Linear solve converged due to CONVERGED_ATOL iterations 1
    4 SNES Function norm 1.461994548486e-03
    Linear solve converged due to CONVERGED_ATOL iterations 1
    5 SNES Function norm 6.066575419934e-04
    Linear solve converged due to CONVERGED_ATOL iterations 1
    6 SNES Function norm 2.655403615568e-04
    Linear solve converged due to CONVERGED_ATOL iterations 1
    7 SNES Function norm 1.202717080189e-04
    Linear solve converged due to CONVERGED_ATOL iterations 1
    8 SNES Function norm 5.567347859811e-05
    Linear solve converged due to CONVERGED_ATOL iterations 1
    9 SNES Function norm 2.613778728633e-05
    Linear solve converged due to CONVERGED_ATOL iterations 1
   10 SNES Function norm 1.238875375025e-05
    Linear solve converged due to CONVERGED_ATOL iterations 1
   11 SNES Function norm 5.911640181274e-06
    Linear solve converged due to CONVERGED_ATOL iterations 1
   12 SNES Function norm 2.834946322421e-06
    Linear solve converged due to CONVERGED_ATOL iterations 1
   13 SNES Function norm 1.364700102970e-06
    Linear solve converged due to CONVERGED_ATOL iterations 1
   14 SNES Function norm 6.589360471852e-07
    Linear solve converged due to CONVERGED_ATOL iterations 1
   15 SNES Function norm 3.189484959867e-07
    Linear solve converged due to CONVERGED_ATOL iterations 1
   16 SNES Function norm 1.547003895237e-07
    Linear solve converged due to CONVERGED_ATOL iterations 1
   17 SNES Function norm 7.516606485367e-08
    Linear solve converged due to CONVERGED_ATOL iterations 0
   18 SNES Function norm 3.657704494472e-08
    Linear solve converged due to CONVERGED_ATOL iterations 0
   19 SNES Function norm 1.782258437327e-08
    Linear solve converged due to CONVERGED_ATOL iterations 1
   20 SNES Function norm 8.694405655955e-09
    Linear solve converged due to CONVERGED_ATOL iterations 0
   21 SNES Function norm 4.245852416264e-09
  Nonlinear solve converged due to CONVERGED_SNORM_RELATIVE iterations 21
1 TS dt 1. time 0.
    0 SNES Function norm 1.826079838384e+01
    Linear solve converged due to CONVERGED_ATOL iterations 31
    1 SNES Function norm 4.299962863814e-02
    Linear solve converged due to CONVERGED_ATOL iterations 1
    2 SNES Function norm 6.431697063491e-03
    Linear solve converged due to CONVERGED_ATOL iterations 1
    3 SNES Function norm 1.639359669638e-03
    Linear solve converged due to CONVERGED_ATOL iterations 1
    4 SNES Function norm 5.026774903429e-04
    Linear solve converged due to CONVERGED_ATOL iterations 1
    5 SNES Function norm 1.680648206351e-04
    Linear solve converged due to CONVERGED_ATOL iterations 1
    6 SNES Function norm 5.967136123324e-05
    Linear solve converged due to CONVERGED_ATOL iterations 1
    7 SNES Function norm 2.230420331383e-05
    Linear solve converged due to CONVERGED_ATOL iterations 1
    8 SNES Function norm 8.721830601750e-06
    Linear solve converged due to CONVERGED_ATOL iterations 1
    9 SNES Function norm 3.542878294948e-06
    Linear solve converged due to CONVERGED_ATOL iterations 1
   10 SNES Function norm 1.484027205445e-06
    Linear solve converged due to CONVERGED_ATOL iterations 1
   11 SNES Function norm 6.368036188603e-07
    Linear solve converged due to CONVERGED_ATOL iterations 1
   12 SNES Function norm 2.784493002713e-07
    Linear solve converged due to CONVERGED_ATOL iterations 1
   13 SNES Function norm 1.235724894148e-07
    Linear solve converged due to CONVERGED_ATOL iterations 1
   14 SNES Function norm 5.549400015107e-08
    Linear solve converged due to CONVERGED_ATOL iterations 0
   15 SNES Function norm 2.516326072475e-08
    Linear solve converged due to CONVERGED_ATOL iterations 0
   16 SNES Function norm 1.150182762543e-08
    Linear solve converged due to CONVERGED_ATOL iterations 1
   17 SNES Function norm 5.293129922422e-09
  Nonlinear solve converged due to CONVERGED_SNORM_RELATIVE iterations 17

# -- many lines omitted --

50 TS dt 1. time 49.
    0 SNES Function norm 1.583816680399e-01
    Linear solve converged due to CONVERGED_ATOL iterations 1
    1 SNES Function norm 3.297785796277e-05
    Linear solve converged due to CONVERGED_ATOL iterations 1
    2 SNES Function norm 1.625267288996e-05
    Linear solve converged due to CONVERGED_ATOL iterations 1
    3 SNES Function norm 8.011050419800e-06
    Linear solve converged due to CONVERGED_ATOL iterations 1
    4 SNES Function norm 3.949253118601e-06
    Linear solve converged due to CONVERGED_ATOL iterations 1
    5 SNES Function norm 1.947152025331e-06
    Linear solve converged due to CONVERGED_ATOL iterations 1
    6 SNES Function norm 9.601580478216e-07
    Linear solve converged due to CONVERGED_ATOL iterations 1
    7 SNES Function norm 4.735242133215e-07
    Linear solve converged due to CONVERGED_ATOL iterations 1
    8 SNES Function norm 2.335591654437e-07
    Linear solve converged due to CONVERGED_ATOL iterations 1
    9 SNES Function norm 1.152140892580e-07
    Linear solve converged due to CONVERGED_ATOL iterations 0
   10 SNES Function norm 5.684167635325e-08
    Linear solve converged due to CONVERGED_ATOL iterations 0
   11 SNES Function norm 2.804655474534e-08
    Linear solve converged due to CONVERGED_ATOL iterations 0
   12 SNES Function norm 1.384019346908e-08
    Linear solve converged due to CONVERGED_ATOL iterations 0
   13 SNES Function norm 6.830511841192e-09
  Nonlinear solve converged due to CONVERGED_SNORM_RELATIVE iterations 13
51 TS dt 1. time 50.
 >> /Users/baagaard/software/unix/py310-venv/pylith-debug/lib/python3.10/site-packages/pylith/problems/Problem.py:204:finalize
 -- timedependent(info)
 -- Finalizing problem.


```

In contrast with Step 1, the evolution of porosity in Step 2 results in nonlinear governing equations.
This requires multiple iterations of the linear solve for the nonlinear solver to converge.

## Visualizing the results

In {numref}`fig:example:magma:2d:step02:solution` we use ParaView to visualize the evolution of the y displacement component using the `viz/plot_dispwarp.py` Python script.
First, we start ParaView from the `examples/magma-2d` directory.

```{code-block} console
---
caption: Open ParaView using the command line.
---
$ PATH_TO_PARAVIEW/paraview

# For macOS, it will be something like
$ /Applications/ParaView-5.10.1.app/Contents/MacOS/paraview
```

Next, we override the default name of the simulation file with the name of the current simulation.

```{code-block} python
---
caption: Set the simulation in the ParaView Python Shell.
---
>>> SIM = "step02_inflation_statevars"
```

Next we run the `viz/plot_dispwarp.py` Python script as described in {ref}`sec-paraview-python-scripts`.

:::{figure-md} fig:example:magma:2d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2 at t=100 yr. The colors indicate the fluid pressure, and the deformation is exaggerated by a factor of 1000." width="75%"/>

Solution for Step 2 at t=100 yr.
The colors of the shaded surface indicate the fluid pressure, and the deformation is exaggerated by a factor of 1000.
The reference state gives rise to greater vertical deformation compared to Step 1.
The choice of material properties does not lead to significant changes in porosity in either material during the simulation (not shown in the figure).
:::
