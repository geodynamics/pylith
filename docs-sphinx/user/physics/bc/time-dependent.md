(sec-user-physics-time-dependent-bc)=
# Time-Dependent Boundary Conditions

Several boundary conditions use a common formulation for the spatial and temporal variation of the boundary condition parameters,
%
\begin{equation}
f(\vec{x})=f_{0}(\vec{x})+\dot{f}_{1}(\vec{x})(t-t_{1}(\vec{x}))+f_{2}(\vec{x})a(t-t_{2}(\vec{x})),\
\end{equation}
%
where

:$f(\vec{x})$: may be a scalar or vector parameter,
:$f_{0}(\vec{x})$: is a constant value (independent of time),
:$\dot{f}_{1}(\vec{x})$: is a constant rate of change in the value with time,
:$t_{1}(\vec{x})$: is the onset time for the constant rate of change,
:$f_{2}(\vec{x})$: is the amplitude for the temporal modulation,
:$a(t)$: is the variation in amplitude with time,
:$t_{2}(\vec{x})$: is the onset time for the temporal modulation, and
:$\vec{x}$: is the position of a location in space.

This common formulation permits easy specification of a scalar or vector with a constant value, constant rate of change of a value, or modulation of a value in time.
One can specify just the initial value, just the rate of change of the value (along with the corresponding onset time), or just the modulation in amplitude (along with the corresponding temporal variation and onset time), or any combination of the three.

## Time-Dependent Dirichlet Boundary Conditions

You can use the `DirichletTimeDependent` boundary conditions to prescribe a solution subfield on a boundary of the finite-element mesh.
The spatial database files for the Dirichlet boundary condition specify the parameters for the time-dependent expression.

:::{important}
The spatial database files for Dirichlet boundary conditions must contain values for all degrees of freedom (x and y for 2-D, and x, y, and z for 3-D even if they are not constrained. This limitation is imposed by the PETSc `DMPlex` interface.
:::

```{table} Values in the auxiliary field spatial databases used for Dirichlet time-dependent boundary conditions.
:name: tab:dirichlet:spatial:database
| Flag               | Required Values                                                     |
|:-------------------|:--------------------------------------------------------------------|
| `use_initial`      | `initial_amplitude_x`, `initial_amplitude_y`, `initial_amplitude_z` |
| `use_rate`         | `rate_start_time`, `rate_amplitude_x`, `rate_amplitude_y`, `rate_amplitude_z` |
| `use_time_history` | `time_history_start`, `time_history_amplitude_x`, `time_history_amplitude_y`, `time_history_amplitude_z` |
```

:::{seealso}
[`DirichletTimeDependent` Component](../../components/bc/DirichletTimeDependent.md) for the Pyre properties and facilities and configuration examples.
:::

## Neumann Time-Dependent Boundary Conditions

Neumann boundary conditions are surface tractions applied over a boundary.
As with the `DirichletTimeDependent` condition, each Neumann boundary condition can only be applied to a simply-connected surface on an external boundary.
The spatial database file the auxiliary subfields for the Neumann boundary condition specify the parameters for the time-dependent expressions.

```{table} Values in the auxiliary field spatial database used for Neumman time-dependent boundary conditions.
| Dimension |    Flag |                       Required Values                                              |
|:----|:----|:----------------------------------------------------------------------------------------------------------------------------|
| 2   | `use_initial`      | `initial_amplitude_normal`,  `initial_amplitude_tangential` |
|     | `use_rate`         | `rate_start_time`, `rate_amplitude_normal`, `rate_amplitude_tangential` |
|     | `use_time_history` | `time_history_start`, `time_history_amplitude_normal`, `time_history_amplitude_tangential` |
| 3   | `use_initial`      | `initial_amplitude_normal`, `initial_amplitude_tangential_1`, `initial_amplitude_tangential_2` |
|     | `use_rate`         | `rate_start_time`, `rate_amplitude_normal`, `rate_amplitude_tangential_1`, `rate_amplitude_tangential_2` |
|     | `use_time_history` | `time_history_start`, `time_history_amplitude_normal`, `time_history_amplitude_tangential_1`, `time_history_amplitude_tangential_2` |
```

:::{seealso}
[`NeumannTimeDependent` Component](../../components/bc/NeumannTimeDependent.md) for the Pyre properties and facilities and configuration examples.
:::
