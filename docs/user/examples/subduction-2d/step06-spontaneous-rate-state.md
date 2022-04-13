# Step 6: Spontaneous Earthquakes With Rate-State Friction

In this example we replace the slip-weakening in Step 5 with rate- and state-friction using the ageing law.
We also lengthen the duration of the simulation to 200 years and reduce the time step to 1.0 years, which were determined through trial and error to get a couple earthquake cycles with reasonable convergence for this relatively coarse resolution mesh.

```{code-block} cfg
---
caption: Excerpt from `step06.cfg`
---

[pylithapp.timedependent.formulation.time_step]
total_time = 200.0*year
dt = 1.0*year
```

The specification of the parameters for the rate- and state-friction model follow a similar pattern to the ones for the slip-weakening friction in Step 5.
Our regularization of the coefficient of friction for near zero slip rate values involves a transition to a linear dependence on slip rate; in this example we specify that this transition should occur at a nondimensional slip rate of $1.0 \times 10^{-6}$.
We impose depth variation of the friction model parameters via a `SimpleGridDB` spatial database in order to generate earthquake-like ruptures in the seismogenic zone with stable sliding above and below.
For the initial tractions, we impose uniform values using a `SimpleDB` spatial database.
We set the initial state for the friction model to be roughly consistent with steady state sliding at the reference coefficient of friction at the reference slip rate, and include it in the state variable in the output as a check.

```{code-block} cfg
---
caption: Excerpt from `step 06.cfg`
---
[pylithapp.timedependent.interfaces.fault_slabtop]
# --- Skipping parameters discussed in previous examples. ---
# Friction
friction = pylith.friction.RateStateAgeing
friction.label = Rate-state friction
# Nondimensional slip rate below which friction depends linearly on slip rate.
friction.linear_slip_rate = 1.0e-6

# Set spatial database for distribution of friction parameters
friction.db_properties = spatialdata.spatialdb.SimpleGridDB
friction.db_properties.label = Slip weakening
friction.db_properties.filename = fault_slabtop_ratestate.spatialdb

# Set spatial database for the initial value of the state variable.
friction.db_initial_state = spatialdata.spatialdb.UniformDB
friction.db_initial_state.label = Rate State Ageing State
friction.db_initial_state.values = [state-variable]
# theta_ss = characteristic_slip_dist / reference_slip_rate
friction.db_initial_state.data = [20.0*year]

# Initial fault tractions
traction_perturbation = pylith.faults.TractPerturbation
traction_perturbation.db_initial = spatialdata.spatialdb.UniformDB
traction_perturbation.db_initial.label = Initial fault tractions
traction_perturbation.db_initial.values = [traction-shear, traction-normal]
traction_perturbation.db_initial.data = [-12.0*MPa, -20.0*MPa]

[pylithapp.problem.interfaces.fault_slabtop.output]
writer = pylith.meshio.DataWriterHDF5
writer.filename = output/step06-fault-slabtop.h5
vertex_info_fields = [normal_dir, strike_dir]
vertex_data_fields = [slip, slip_rate, traction, state_variable]
```

```{code-block} console
---
caption: Run Step 6 simulation
---
$ pylith step06.cfg
```

The problem will produce fourteen pairs of HDF5/Xdmf files.
{numref}`fig:example:subduction:2d:step06`, which was created using the ParaView Python script `viz/plot_dispwarp.py`, displays the magnitude of the velocity field with the original configuration exaggerated by a factor of 4000.
Steady slip is largely confined to the stable sliding regions with a sequence of ruptures in the seismogenic zone; note how the rate-state friction allows a more natural nucleation of the ruptures compared to the slip-weakening friction.
{numref}`fig:example:subduction:2d:step06:slip` shows the cumulative slip as a function of time and distance down dip from the trench.

:::{figure-md} fig:example:subduction:2d:step06
<img src="figs/step06_soln.*" alt="Solution for Step 6 at the end of the simulation. The colors indicate the magnitude of the x-displacement component and the deformation has been exaggerated by a factor of 10,000." width="100%"/>

Solution for Step 6 at the end of the simulation. The colors indicate the magnitude of the x-displacement component and the deformation has been exaggerated by a factor of 10,000.
:::

:::{figure-md} fig:example:subduction:2d:step06:slip
<img src="figs/step06_slip.*" alt="Cumulative slip as a function of time and depth in Step 6. The red lines indicate slip every 10 time steps." width = "100%"/>

Cumulative slip as a function of time and depth in Step 6. The red lines indicate slip every 10 time steps.
:::
