# Step 7: Slip on Two Faults and Maxwell Viscoelastic Materials

% Metadata extracted from parameter files.
```{include} step07_twofaults_maxwell-synopsis.md
```

## Simulation parameters

In this example we replace the linear, isotropic elastic bulk rheology in the slab with a linear, isotropic Maxwell viscoelastic rheology.
We use the same boundary conditions as in Step 6 but reduce the time step to resolve the viscoelastic relaxation.
The parameters specific to this example are in `step07_twofaults_maxwell.cfg`.

We use a very short relaxation time of 20 years, so we run the simulation for 100 years with a time step of 4 years.
We use a starting time of -4 years so that the first time step will advance the solution time to 0 years.

```{code-block} cfg
---
caption: Time stepping parameters for Step 7.
---
[pylithapp.problem]
initial_dt = 4.0*year
start_time = -4.0*year
end_time = 100.0*year

scales.relaxation_time = 20.0*year
```

The uniform slip creates a large stress concentration at the end of the fault.
We improve the resolution by using uniform global refinement and discretizing the displacement and Lagrange multiplier fields with a basis order of 2.

```{code-block} cfg
---
caption: Enabling uniform global refinement and discretizing the displacement and Lagrange multiplier fields with a basis order of 2 in Step 7.
---
[pylithapp.mesh_generator]
refiner = pylith.topology.RefineUniform

[pylithapp.problem]
solution = pylith.problems.SolnDispLagrange
defaults.quadrature_order = 2

[pylithapp.problem.solution.subfields]
displacement.basis_order = 2
lagrange_multiplier_fault.basis_order = 2
```

```{code-block} cfg
---
caption: Maxwell viscoelastic bulk rheology parameters for the slab in Step 7.
---
[pylithapp.problem.materials]
slab.bulk_rheology = pylith.materials.IsotropicLinearMaxwell

[pylithapp.problem.materials.slab]
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Maxwell viscoelastic properties
db_auxiliary_field.iohandler.filename = mat_maxwell.spatialdb

bulk_rheology.auxiliary_subfields.maxwell_time.basis_order = 0
```

## Running the simulation

```{code-block} console
---
caption: Run Step 7 simulation
---
$ pylith step07_twofaults_maxwell.cfg

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
    0 SNES Function norm 1.841529893626e-01
      Linear solve converged due to CONVERGED_ATOL iterations 16
    1 SNES Function norm 4.357719233336e-09
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
26 TS dt 0.2 time 5.
 >> /software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:232:finalize
 -- timedependent(info)
 -- Finalizing problem.```

From the end of the output written to the terminal window, we see that the simulation advanced the solution 26 time steps.
The PETSc TS display time in the nondimensional units, so a time of 5 corresponds to 100 years.

## Visualizing the results

In {numref}`fig:example:reverse:2d:step07:solution` we use the `pylith_viz` utility to visualize the x displacement field.
You can move the slider or use the `p` and `n` keys to change the increment or decrement time.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filename=output/step07_twofaults_maxwell-domain.h5 warp_grid --component=x
```

:::{figure-md} fig:example:reverse:2d:step07:solution
<img src="figs/step07-solution.*" alt="Solution for Step 7 at t=100 yr. The colors indicate the x displacement, and the deformation is exaggerated by a factor of 1000." width="600px"/>

Solution for Step 7 at t=100 years.
The colors of the shaded surface indicate the x displacement, and the deformation is exaggerated by a factor of 1000.
The undeformed configuration is shown by the gray wireframe.
Viscoelastic relaxation results in significant deformation in the `slab` material.
:::

## Step 7b: Adaptive Time Stepping

In Step 7b we demonstrate the use adaptive time stepping.
We start with an initial time step of 0.2 years and let the adaptive time stepping algorithm increase the time step as the rate of deformation decreases.

```{code-block} console
---
caption: Run Step 7b simulation
---
$ pylith step07_twofaults_maxwell.cfg step07b_twofaults_maxwell.cfg

# The output should look something like the following.
# -- many lines omitted --

30 TS dt 0.120677 time 4.87932
    0 SNES Function norm 1.161325530040e-01
      Linear solve converged due to CONVERGED_ATOL iterations 19
    1 SNES Function norm 8.829574302098e-09
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
      TSAdapt basic beuler 0: step  30 accepted t=4.87932    + 1.207e-01 dt=2.739e-01  wlte=0.00776  wltea=   -1 wlter=   -1
31 TS dt 0.273922 time 5.
 >> /Users/brad/software/unix/py3.12-venv/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:232:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

With smaller time steps following the two earthquakes, this simulation takes 5 more time steps than Step 7 with uniform time steps.
For simulations covering longer periods of time or a single earthquake, the adaptive time stepping usually results in substantially fewer time steps.

```{code-block} console
---
caption: Compare the results from Step 7 ad 7b using the `plot_compare.py` Python script.
---
./plot_compare.py --step=7
```

:::{figure-md} fig:example:reverse:2d:step07b:solution
<img src="figs/step07-compare.*" alt="Displacement and shear stress time histories for Step 7 and 7b." width="600px"/>

Displacement and shear stress time histories for Step 7 and 7b at a location in the viscoelastic slab below the bottom ends of the main fault and splay fault.
The backward Euler time stepping algorithm does a poor job of resolving the rapid deformation following the sudden application of coseismic slip, resulting in some oscillations in the deformation in Step 7b which uses adaptive time stepping.
:::
