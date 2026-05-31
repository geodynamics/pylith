# Step 2: Gravitational Body Forces with Reference Stress

% Metadata extracted from parameter files.
```{include} step02_gravity_refstate-synopsis.md
```

## Simulation parameters

This example involves using a reference stress state to minimize the deformation when we apply the gravitational body forces.
The solution will be the perturbation from the reference state with zero displacements.
This is one method for obtaining an initial stress state associated with gravitational body forces.
We use the same roller boundary conditions that we used in Step 1.
The parameters specific to this example are in `step02_gravity_refstate.cfg`.

We use a reference stress state that matches the overburden (lithostatic) pressure.
We have uniform material properties, so the overburden is
%
\begin{equation}
\sigma_{xx} = \sigma_{yy} = \sigma_{zz} = \int_0^z \rho g \, dz = \rho g z,
\end{equation}
where compressive stress is negative.

```{code-block} cfg
---
caption: Parameters for reference stresses for Step 2. We only show the details for the slab material.
---
[pylithapp.problem.materials.slab]
db_auxiliary_field.iohandler.filename = mat_gravity_refstate.spatialdb
db_auxiliary_field.query_type = linear

[pylithapp.problem.materials.slab.bulk_rheology]
use_reference_state = True

auxiliary_subfields.reference_stress.basis_order = 1
auxiliary_subfields.reference_strain.basis_order = 0
```

## Running the simulation

```{code-block} console
---
caption: Run Step 2 simulation
---
$ pylith step02_gravity_refstate.cfg

# The output should look something like the following.
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/apps/PyLithApp.py:79:main
 -- info (application-flow)
 -- Running on 1 process(es).

# -- many lines omitted --

 >> src/cig/pylith/libsrc/pylith/problems/TimeDependent.cc:473:void pylith::problems::TimeDependent::solve()
 -- info (application-flow)
 -- Component 'timedependent.problem': Solving equations.
0 TS dt 0.001 time 0.
    0 SNES Function norm 6.916025621683e-12
    Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 0
1 TS dt 0.001 time 0.001
 >> software/pylith-debug/lib/python3.12/site-packages/pylith/problems/Problem.py:222:finalize
 -- info (application-flow)
 -- Finalizing problem.
WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
There is one unused database option. It is:
Option left: name:-mg_fine_ksp_max_it value: 5 source: code
```

By design we set the reference stress state so that it matches the loading from gravitational body forces in our domain with uniform material properties.
As a result, the first nonlinear solver residual evaluation meets the convergence criteria.
The linear solver is not used; this is why PETSc reports an unused option at the end of the simulation.

## Visualizing the results

In {numref}`fig:example:reverse:2d:step02:solution` and {numref}`fig:example:reverse:2d:step02:stress` we use the `pylith_viz` utility to visualize the simulation results.

```{code-block} console
---
caption: Visualize PyLith output using `pylith_viz`.
---
pylith_viz --filenames=output/step02_gravity_refstate-domain.h5 warp_grid --exaggeration=5
pylith_viz --filenames=output/step02_gravity_refstate-crust.h5,output/step02_gravity_refstate-slab.h5,output/step02_gravity_refstate-wedge.h5 warp_grid --field=cauchy_stress --component=xy --exaggeration=5
```

:::{figure-md} fig:example:reverse:2d:step02:solution
<img src="figs/step02-solution.*" alt="Solution for Step 2. The colors indicate the magnitude of the displacement." width="600px"/>

Solution for Step 2.
The colors of the shaded surface indicate the magnitude of the displacement, which is zero.
The undeformed configuration is shown by the gray wireframe.
:::

:::{figure-md} fig:example:reverse:2d:step02:stress
<img src="figs/step02-stress.*" alt="Cauchy stress tensor component xy for Step 2. The colors indicate the stress tensor component, and the deformation is exaggerated by a factor of 5." width="600px"/>

Cauchy stress tensor component xy for Step 2.
The colors of the shaded surface indicate the xy component of the Cauchy stress tensor, and the deformation is exaggerated by a factor of 5.
The undeformed configuration is shown by the gray wireframe.
The shear stress is zero.
:::
