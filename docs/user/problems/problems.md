# Types of Simulations

PyLith currently supports two types of problems:

* Time dependent, and
* Static Green's functions.

(sec-user-problems-time-dependent)=
## Time-Dependent Problem (`TimeDependent`)

This type of problem solves a time-dependent boundary value problem and applies to static, quasistatic, and dynamic simulations.
We use the PETSc object to manage the time-stepping, including the selection of the time-stepping algorithm.
By default PyLith uses the backward Euler time-stepping algorithm.

:::{admonition} Pyre User Interface
:class: seealso
See [`TimeDependent` Component](../components/problems/TimeDependent.md) for Pyre properties and facilities and configuration examples.
:::

### Initial Conditions

:::{note}
New in v3.0.0.
:::

The initial conditions for a simulation are specified via a combination of initial values for the solution and initial values for state variables.
The initial values for state variables are specified via the spatial databases for the auxiliary field of each material.
In this section we discuss how to set the initial values of the solution field.

#### `InitialConditionDomain`

We use this object when we want to specify the initial values of solution subfields across the entire domain using a single spatial database.

:::{admonition} Pyre User Interface
:class: seealso
See [`InitialConditionDomain` Component](../components/problems/InitialConditionDomain.md) for Pyre properties and facilities and configuration examples.
:::

#### `InitialConditionPatch`

We use this object when we want to specify the initial values of solution subfields across patches of the domain defined by materials.
Each patch specifies the initial values for a single material.
Not all materials need to have initial values.

:::{admonition} Pyre User Interface
:class: seealso
See [`InitialConditionPatch` Component](../components/problems/InitialConditionPatch.md) for Pyre properties and facilities and configuration examples.
:::

### Numerical Damping in Explicit Time Stepping

:::{danger}
Not yet reimplemented in v3.x.
:::

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

(sec-user-problems-greensfns)=
## Green&rsquo;s Functions Problem (`GreensFns`)

This type of problem applies to computing static Green's functions for elastic deformation.
The `GreensFns` problem loops over a suite of fault slip impulses and computes the static solution for each impulse using the linear solver.
In the output files, the deformation at each "time step" is the deformation for a different slip impulse.
The fault slip impulses are specified using `FaultCohesiveImpulses` for the fault.
See {ref}`sec-user-physics-fault-cohesive-impulses` for more information.

:::{warning}
The `GreensFns` problem generates slip impulses on a fault.
PyLith currently requires that impulses be applied to a single fault of type `FaultCohesiveImpulses`.
:::

:::{admonition} Pyre User Interface
:class: seealso
See [`GreensFns` Component](../components/problems/GreensFns.md) for Pyre properties and facilities and configuration examples.
:::
