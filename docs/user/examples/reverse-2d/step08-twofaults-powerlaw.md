(sec-user-examples-reverse-2d-step08)=
# Step 8: Slip on Two Faults and Power-law Viscoelastic Materials

% Metadata extracted from parameter files.
```{include} step08_twofaults_powerlaw-synopsis.md
```

## Simulation parameters

In this example we replace the linear, isotropic Maxwell viscoelastic bulk rheology in the slab in Step 7 with an isotropic powerlaw viscoelastic rheology.
The other parameters remain the same as those in Step 7.
The parameters specific to this example are in `step08_twofaults_powerlaw.cfg`.

```{code-block} cfg
---
caption: Power-law viscoelastic bulk rheology parameters for Step 8.
---
[pylithapp.problem.materials]
slab.bulk_rheology = pylith.materials.IsotropicPowerLaw

[pylithapp.problem.materials.slab]
db_auxiliary_field = spatialdata.spatialdb.CompositeDB
db_auxiliary_field.description = Power law material properties

bulk_rheology.auxiliary_subfields.power_law_reference_strain_rate.basis_order = 0
bulk_rheology.auxiliary_subfields.power_law_reference_stress.basis_order = 0
bulk_rheology.auxiliary_subfields.power_law_exponent.basis_order = 0

[pylithapp.problem.materials.slab.db_auxiliary_field]
# Elastic properties
values_A = [density, vs, vp]
db_A = spatialdata.spatialdb.SimpleDB
db_A.description = Elastic properties for slab
db_A.iohandler.filename = mat_elastic.spatialdb

# Power law properties
values_B = [
	 power_law_reference_stress,  power_law_reference_strain_rate,  power_law_exponent,
	 viscous_strain_xx, viscous_strain_yy, viscous_strain_zz, viscous_strain_xy,
	 reference_stress_xx, reference_stress_yy, reference_stress_zz, reference_stress_xy,
	 reference_strain_xx, reference_strain_yy, reference_strain_zz, reference_strain_xy,
	 deviatoric_stress_xx,  deviatoric_stress_yy,  deviatoric_stress_zz,  deviatoric_stress_xy
	 ]
db_B = spatialdata.spatialdb.SimpleDB
db_B.description = Material properties specific to power law bulk rheology for the slab
db_B.iohandler.filename = mat_powerlaw.spatialdb
db_B.query_type = linear
```

## Power-law spatial database

*New in v4.0.0.*

We use the utility script `pylith_powerlaw_gendb` (see {ref}`sec-user-run-pylith-pylith-powerlaw-gendb`) to generate the spatial database `mat_powerlaw.spatialdb` with the power-law bulk rheology parameters.
We provide the parameters for `pylith_powerlaw_gendb` in `powerlaw_gendb.cfg`, which follows the same formatting conventions as the PyLith parameter files. 

```{code-block} console
---
caption: Generate spatial database with power-law viscoelastic material properties.
---
pylith_powerlaw_gendb powerlaw_gendb.cfg
```

## Running the simulation

```{code-block} console
---
caption: Run Step 8 simulation
---
$ pylith step08_twofaults_powerlaw.cfg

# The output should look something like the following.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:77:main
 -- pylithapp(info)
 -- Running on 1 process(es).
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/meshio/MeshIOObj.py:38:read
 -- meshiopetsc(info)
 -- Reading finite-element mesh
 >> /src/cig/pylith/libsrc/pylith/meshio/MeshIO.cc:85:void pylith::meshio::MeshIO::read(pylith::topology::Mesh *, const bool)
 -- meshiopetsc(info)
 -- Component 'reader': Domain bounding box:
    (-100000, 100000)
    (-100000, 0)

# -- many lines omitted --

25 TS dt 0.2 time 4.8
    0 SNES Function norm 3.244406041330e-02
      Linear solve converged due to CONVERGED_ATOL iterations 16
    1 SNES Function norm 4.276521889496e-04
      Linear solve converged due to CONVERGED_ATOL iterations 10
    2 SNES Function norm 1.467733852722e-05
      Linear solve converged due to CONVERGED_ATOL iterations 6
    3 SNES Function norm 5.094245096618e-07
      Linear solve converged due to CONVERGED_ATOL iterations 2
    4 SNES Function norm 1.675562142387e-08
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 4
26 TS dt 0.2 time 5.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:232:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

As in Step 7, the simulation advances 26 time steps.
With a nonlinear bulk rheology, the nonlinear solver now requires a few iterations to converge at each time step.

## Visualizing the results

In {numref}`fig:example:reverse:2d:step08:solution` we use the `pylith_viz` utility to visualize the x displacement field.
You can move the slider or use the `p` and `n` keys to change the increment or decrement time.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step08_twofaults_powerlaw-domain.h5 warp_grid --component=x
```

:::{figure-md} fig:example:reverse:2d:step08:solution
<img src="figs/step08-solution.*" alt="Solution for Step 8. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 1000." width="600px"/>

Solution for Step 8 at t=100 years.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is shown by the gray wireframe.
Our parameters for the power-law bulk rheology result in much less viscoelastic relaxation in this case compared to Step 7.
:::

## Step 8b: Adaptive Time Stepping

In Step 8b we demonstrate the use adaptive time stepping.
We start with an initial time step of 0.2 years and let the adaptive time stepping algorithm increase the time step as the rate of deformation decreases.

```{code-block} console
---
caption: Run Step 8b simulation
---
$ pylith step08_twofaults_maxwell.cfg step08b_twofaults_maxwell.cfg

# The output should look something like the following.
# -- many lines omitted --

14 TS dt 1.14678 time 3.85322
    0 SNES Function norm 2.159008163685e-01
      Linear solve converged due to CONVERGED_ATOL iterations 18
    1 SNES Function norm 1.088600033252e-02
      Linear solve converged due to CONVERGED_ATOL iterations 14
    2 SNES Function norm 7.037994748778e-04
      Linear solve converged due to CONVERGED_ATOL iterations 10
    3 SNES Function norm 4.539417964138e-05
      Linear solve converged due to CONVERGED_ATOL iterations 7
    4 SNES Function norm 2.925595088171e-06
      Linear solve converged due to CONVERGED_ATOL iterations 3
    5 SNES Function norm 1.884520992392e-07
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 5
      TSAdapt basic beuler 0: step  14 accepted t=3.85322    + 1.147e+00 dt=2.660e+00  wlte=0.00744  wltea=   -1 wlter=   -1
15 TS dt 2.65981 time 5.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:232:finalize
 -- timedependent(info)
 -- Finalizing problem.```

The deformation rate is smoother with the power-law rheology yielding more efficient results with the adaptive time stepper which uses 11 few time steps compared with the uniform time step of 4 years.

```{code-block} console
---
caption: Compare the results from Step 8 ad 8b using the `plot_compare.py` Python script.
---
./plot_compare.py --step=8
```

:::{figure-md} fig:example:reverse:2d:step08b:solution
<img src="figs/step08-compare.*" alt="Displacement and shear stress time histories for Step 8 and 8b." width="600px"/>

Displacement and shear stress time histories for Step 8 and 8b at a location in the viscoelastic slab below the bottom ends of the main fault and splay fault.
The adaptive time stepping does a better job of capturing the variable rate of deformation associated with the power-law viscoelastic rheology; however, the simulation steps over the time of the earthquake on the splay fault before reducing the time step.
:::
