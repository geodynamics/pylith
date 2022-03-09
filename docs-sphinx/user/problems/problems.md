# Simulation Types (**problem**)

## Time-Dependent Problem (*TimeDependent*)

This type of problem applies to transient static, quasistatic, and dynamic simulations.
We use the PETSc object to manage the time-stepping.
The selection of the time-stepping algorithm is done via the PETSc options.
By default PyLith use the backward Euler time stepping algorithm.

```{code-block} cfg
---
caption: Setting PETSc time-stepping options in a `cfg` file
---
[pylithapp.petsc]
ts_type = beuler
ts_monitor = true
ts_error_if_step_fails = true
```

The *TimeDependent* properties and facilities include

:initial_dt: Initial time step (default=1.0*year);
:start_time: Starting time of the problem (default=0.0\*year, the first time step will be from **start_time** to **start_time** + **initial_step**)

:total_time: Time duration of the problem (default=0.0*year);

:max_timesteps: Maximum number of time steps (default=20000);

:ic: Initial conditions for solution (default=*EmptyBin*); and

:notify_observers_ic: Send observers solution with initial conditions before time stepping (default=False);

```{code-block} cfg
---
caption: *TimeDependent* parameters in a `cfg` file
---
[pylithapp.timedependent]
initial_dt = 10.0*year
start_time = 0.0*year
end_time = 100.0*year
ic = [pylith.problems.InitialConditionDomain]
```

### Initial Conditions

The initial conditions for a simulation are specified via a combination of the initial values for the solution and the initial values for the state variables.
The initial values for state variables are specified via the spatial databases for the auxiliary field.
In this section we discuss how to set the initial values of the solution field.

#### *InitialConditionDomain*

We use this object when we want to specify the initial values of solution subfields across the entire domain using a single spatial database.
The properties and facilities of the are:

:subfields: List of names of solution subfields for which initial values are supplied (default=[displacement]); and
:db: Spatial database with initial values for solution subfields (default=*spatialdata.spatialdb.SimpleDB*)

#### *InitialConditionPatch*

We use this object when we want to specify the initial values of solution subfields across patches of the domain defined by materials.
Each patch specifies the initial values for a single material.
Not all materials need to have initial values.
The properties and facilities of the are:

:subfields: List of names of solution subfields for which initial values are supplied (default=[displacement]); and
:id: Id of material associated with the patch (default=0); and
:db: Spatial database with initial values for solution subfields (default=*spatialdata.spatialdb.SimpleDB*)

## Numerical Damping in Explicit Time Stepping

Not yet reimplemented in v3.x.

%In explicit time-stepping formulations for elasticity, boundary conditions and fault slip can excite short waveform elastic waves that are not accurately resolved by the discretization.
%We use numerical damping via an artificial viscosity {cite}`Knopoff:Ni:2001,Day:Ely:2002` to reduce these high frequency oscillations.
%In computing the strains for the elasticity term in equation \ref{eq:elasticity:integral:dynamic:t}, we use an adjusted displacement rather than the actual displacement, where
%
%```{math}
%\vec{u}^{adj}(t)=\vec{u}(t)+\eta^{*}\Delta t\vec{\dot{u}}(t),
%```
%
%$\vec{u}^{adj}(t)$ is the adjusted displacement at time $t$, $\vec{u}(t)$ is the original displacement at time (t), $\eta^{*}$is the normalized artificial viscosity, $\Delta t$ is the time step, and $\vec{\dot{u}}(t)$ is the velocity at time $t$.
%The default value for the normalized artificial viscosity is 0.1.
%We have found values in the range 0.1-0.4 sufficiently suppress numerical noise while not excessively reducing the peak velocity.
%An example of setting the normalized artificial viscosity in a `cfg` file is
%
%```{code-block} cfg
%[pylithapp.timedependent.formulation]
%norm_viscosity = 0.2
%```
%
## Green&rsquo;s Functions Problem (*GreensFns*)

Not yet reimplemented in v3.x.

%This type of problem applies to computing static Green's functions for elastic deformation.
%The *GreensFns* problem specializes the time-dependent facility to the case of static simulations with slip impulses on a fault.
%The default formulation is the Implicit formulation and should not be changed as the other formulations are not applicable to static Green's functions.
%In the output files, the deformation at each "time step" is the deformation for a different slip impulse.
%The properties provide the ability to select which fault to use for slip impulses.
%The only fault component available for use with the *GreensFns* problem is the *FaultCohesiveImpulses* component discussed in \ref{sec:fault:cohesive:impulses}.
%The *GreensFns* properties amd facilities include:
%
%:fault_id: Id of fault on which to impose slip impulses.
%:formulation: Formulation for solving the partial differential equation.
%:progress_monitor: Simple progress monitor via text file.
%
%```{code-block} cfg
%---
%caption: *GreensFns* parameters in a `cfg` file]
%---
%[pylithapp]
%problem = pylith.problems.GreensFns ; Change problem type from the default
%
%[pylithapp.greensfns]
%fault_id = 100 ; Default value
%formulation = pylith.problems.Implicit ; default
%progres_monitor = pylith.problems.ProgressMonitorTime ; default
%```
%
%:::{warning}
%The *GreensFns* problem generates slip impulses on a fault.
%The current version of PyLith requires that impulses can only be applied to a single fault and the fault facility must be set to *FaultCohesiveImpulses*.
%:::
%
### Progress Monitors

Not yet reimplemented in v3.x.

%The progress monitors make it easy to monitor the general progress of long simulations, especially on clusters where stdout is not always easily accessible.
%The progress monitors update a simulation's current progress by writing information to a text file.
%The information includes time stamps, percent completed, and an estimate of when the simulation will finish.
%
%#### *ProgressMonitorTime*
%
%This is the default progress monitor for time-stepping problems.
%The monitor calculates the percent completed based on the time at the current time step and the total simulated time of the simulation, not the total number of time steps (which may be unknown in simulations with adaptive time stepping).
%The *ProgressMonitorTime} properties include:
%
%:update_percent: Frequency (in percent) of progress updates.
%:filename: Name of output file.
%:t_units: Units for simulation time in output
%
%```{code-block} cfg
%---
%caption: *ProgressMonitorTime* parameters in a `cfg` file
%---
%[pylithapp.problem.progressmonitor]
%update_percent = 5.0 ; default
%filename = progress.txt ; default
%t_units = year ; default
%```
%
%#### *ProgressMonitorStep*
%
%This is the default progress monitor for problems with a specified number of steps, such as Green's function problems.
%The monitor calculates the percent completed based on the number of steps (e.g., Green's function impulses completed).
%The ProgressMonitorStep propertiles include:
%
%:update_percent: Frequency (in percent) of progress updates.
%:filename: Name of output file.
%
%```{code-block} cfg
%---
%caption: *ProgressMonitorSte* parameters in a `cfg` file
%---
%[pylithapp.problem.progressmonitor]
%update_percent = 5.0 ; default
%filename = progress.txt ; default
%```
