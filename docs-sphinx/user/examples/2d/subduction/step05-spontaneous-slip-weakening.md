# Step 5: Spontaneous Earthquakes With Slip-Weakening Friction

We simulate earthquake cycles over 100 years with spontaneous rupture using slip-weakening friction.
As in Step 4 including fault friction requires the nonlinear solver.
Through trial and error we choose a time step of 2.5 years that permits reasonable convergence of the nonlinear solver and runtime.

```{code-block} cfg
---
caption: Excerpt from `step05.cfg`
---
[pylithapp.problem.formulation]
# Fault friction is a nonlinear problem so we need to use the
# nonlinear solver.
solver = pylith.problems.SolverNonlinear

[pylithapp.timedependent.formulation.time_step]
total_time = 100.0*year
dt = 2.5*year
```

In simulations for research purposes, we would use a higher resolution mesh and smaller time steps and investigate the robustness of the solution to these parameters.

We constrain the displacement normal to the lateral and bottom boundaries without restraining the subducting slab.
We also constrain the vertical deformation of the west boundary to facilitate the downward motion of the subducting slab.

```{code-block} cfg
---
caption: Excerpt from `step05.cfg`
---
[pylithapp.timedependent.bc.boundary_west]
bc_dof = [0, 1]
label = bndry_west
db_initial.label = Dirichlet BC on west boundary
```

We replace the prescribed aseismic slip on the subduction interface that we used in Step 2 with a friction interface with the slip-weakening fault constitutive model.

```{code-block} cfg
---
caption: Excerpt from `step05.cfg`
---
[pylithapp.timedependent]
interfaces = [fault_slabtop, fault_slabbot]

# Set the type of fault interface conditions.
[pylithapp.timedependent.interfaces]
fault_slabtop = pylith.faults.FaultCohesiveDyn
fault_slabbot = pylith.faults.FaultCohesiveKin
```

In order to generate stick-slip events, we need the coefficient of friction to decrease with slip.
We choose a slip-weakening friction model with a dynamic coefficient of friction that is less than the static coefficient of friction to provide this behavior.
In quasistatic modeling we use time steps much longer than the slip rise time in an earthquake, so we want the slip confined to one time step or just a few time steps.
This means the drop in the coefficient of friction should be independent in each time step; that is, we want the fault to fully heal between time steps.
This corresponds to setting the **force_healing** property of the *SlipWeakening* object.

A common feature in numerical modeling of subduction zones is stable sliding near the trench and below the seismogenic zone.
We implement stable sliding with the slip-weakening friction via a constant coefficient of friction (equal values for the static and dynamic coefficients of friction).
We create a lower dynamic coefficient of friction in the seismogenic zone, by introducing depth-dependent variations in the dynamic coefficient of friction using a *SimpleGridDB* spatial database as discussed in {ref}`sec:spatial:databases`.
This provides more efficient interpolation compared to the *SimpleDB* implementation.
We impose initial tractions on the fault in a similar fashion as we did in Step 4.
We reduce the initial shear tractions slightly in the seismogenic zone, consistent with a stress drop in the penultimate earthquake followed by loading during the interseismic period.

```{code-block} cfg
---
caption: Excerpt from `step05.cfg`
---
[pylithapp.timedependent.interfaces.fault_slabtop]
# --- Skipping general information discussed previously ---
# Friction
friction = pylith.friction.SlipWeakening
friction.label = Slip weakening
# Force healing after each time step, so weakening is confined to each
# time step and is not carried over into subsequent time steps.
friction.force_healing = True

friction.db_properties = spatialdata.spatialdb.SimpleGridDB
friction.db_properties.label = Slip weakening
friction.db_properties.filename = fault_slabtop_slipweakening.spatialdb

# Initial fault tractions
traction_perturbation = pylith.faults.TractPerturbation
traction_perturbation.db_initial = spatialdata.spatialdb.SimpleGridDB
traction_perturbation.db_initial.label = Initial fault tractions
traction_perturbation.db_initial.filename = fault_slabtop_tractions.spatialdb
```

We adjust several of the solver tolerances.
In general, we impose larger tolerances to reduce runtime at the expense of a less accurate solution.
We set the zero tolerances for detecting slip and suppressing fault opening to {math}`1.0\times 10^{-8}`.
We want tolerances for the linear solve to be smaller than these values, so we use an absolute tolerance of {math}`1.0\times 10^{-9}` and a very small relative tolerance to force the residual below the absolute tolerance.
We impose an absolute tolerance for the nonlinear solver to be greater than our zero tolerances and also force the residual to match the absolute tolerance level by using a very small relative tolerances.
Finally, we set the parameters for the solver used to calculate consistent values for the change in slip for a given change in the Lagrange multipliers (which we sometimes call the friction sensitivity solve).

```{code-block} cfg
---
caption: Excerpt from `step05.cfg`
---
[pylithapp.timedependent.interfaces.fault_slabtop]
zero_tolerance = 1.0e-8
zero_tolerance_normal = 1.0e-8

# Convergence parameters.
ksp_rtol = 1.0e-20
ksp_atol = 1.0e-9
ksp_max_it = 1000

snes_rtol = 1.0e-20
snes_atol = 1.0e-7
snes_max_it = 1000

# Friction sensitivity solve used to compute the increment in slip
# associated with changes in the Lagrange multiplier imposed by the
# fault constitutive model.
friction_pc_type = asm
friction_sub_pc_factor_shift_type = nonzero
friction_ksp_max_it = 25
friction_ksp_gmres_restart = 30
friction_ksp_error_if_not_converged = true
```

```{code-block} console
---
caption: Run Step 5 simulation
---
$ pylith step05.cfg
```

The problem will produce fourteen pairs of HDF5/Xdmf files.
{numref}`fig:example:subduction:2d:step05`, which was created using the ParaView Python script `viz/plot_dispwarp.py` (see {ref}`sec:ParaView:Python:scripts` for a discussion of how to run ParaView Python scripts), displays the magnitude of the velocity field with the original configuration exaggerated by a factor of 4000.
Steady slip is largely confined to the stable sliding regions with a sequence of ruptures in the seismogenic zone; most have a duration of a few time steps, although most of the slip occurs in a single time step.
{numref}`fig:example:subduction:2d:step05:slip` shows the cumulative slip as a function of time and distance down dip from the trench.

:::{figure-md} fig:example:subduction:2d:step05
<img src="figs/step05_soln.*" alt="Solution for Step 5 at the end of the simulation. The colors indicate the magnitude of the x-displacement component and the deformation has been exaggerated by a factor of 10,000. " width="100%" />

Solution for Step 5 at the end of the simulation. The colors indicate the magnitude of the x-displacement component and the deformation has been exaggerated by a factor of 10,000.
:::

:::{figure-md} fig:example:subduction:2d:step05:slip
<img src="figs/step05_slip.*" alt="Cumulative slip as a function of time and depth in Step 5. The red lines indicate slip every 10 time steps." width = "100%" />

Cumulative slip as a function of time and depth in Step 5. The red lines indicate slip every 10 time steps.
:::
