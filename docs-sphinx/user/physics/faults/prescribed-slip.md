(sec-user-physics-prescribed-slip)=
# Prescribed Slip (`FaultCohesiveKin`)

Prescribed slip, often called kinematic earthquake ruptures, use the `FaultCohesiveKin` component to prescribe the slip as a function of time on the fault surface.
Slip may evolve simultaneously over the fault surface instantaneously in a single time step (as is usually done in quasistatic simulations) or propagate over the fault surface over hundreds and up to thousands of time steps (as is usually done in a dynamic simulation).

Multiple earthquake ruptures can be specified on a single fault surface.
This permits repeatedly rupturing the same portion of a fault or combining earthquake rupture on one subset of the fault surface with steady aseismic slip on another subset (the two subsets may overlap in both time and space).
An array of kinematic earthquake rupture components associates a name (string) with each kinematic rupture.
The default dynamic array contains a single earthquake rupture, `rupture`.

:::{admonition} TODO
:class: error

Add note about discretization of "slip" auxiliary subfield.
:::

## Prescribed Slip Parameters (`KinSrc`)

The kinematic rupture parameters include the origin time and slip time function.
The slip initiation time in the slip time function is relative to the origin time (default is 0).
This means that slip initiates at a point at a time corresponding to the sum of the kinematic rupture's origin time and the slip initiation time for that point.

PyLith supports specification of the evolution of fault slip using analytical expressions for the slip time history at each point, where the parameters for the slip time function may vary over the fault surface.
The following slip time functions are available:

:`KinSrcStep`: a step function for quasistatic modeling of earthquake rupture,
:`KinSrcConstRate`: a constant rate function for steady slip (for example, creep),
:`KinSrcRamp`: a constant rate function for a specified time interval,
:`KinSrcBrune`: a slip time function corresponding to the integral of Brune's far-field time function for dynamic modeling of earthquake rupture,
:`KinSrcLiuCos`: a slip time function built on cosine functions for dynamic modeling of earthquake rupture, and
:`KinSrcTimeHistory`: a user-specified slip time function.

### Step-Function Slip Time Function (`KinSrcStep`)

This slip function prescribes a step in slip at a given time at a point:
%
\begin{gather}
D(t)=\left\{ \begin{array}{cc}
0 & 0\leq t<t_0\\
D_{final} & t\ge t_0
\end{array}\right.\,,
\end{gather}
%
where $D(t)$ is slip at time $t$, $D_{final}$ is the final slip, and $t_0$ is the slip initiation time (time when rupture reaches the location).
The slip is specified independently for each of the components of slip, and the slip and slip starting time may vary over the fault surface.

```{table} Values in the auxiliary field spatial database for KinSrcStep.
:name: tab:slip:function:step
|      Subfield     |           Components                 |
|:------------------|:-------------------------------------|
| `initiation_time` |     --                               |
| `final_slip`      | `opening`, `left_lateral`, `reverse` |
```

### Constant Slip Rate Slip Time Function (`KinSrcConstRate`)

This slip function prescribes a constant slip rate for the evolution of slip at a point:
%
\begin{gather}
  D(t)=\left\{ \begin{array}{cc}
0 & 0\leq t<t_0\\
V(t-t_0) & t\ge t_0
\end{array}\right.\,,
\end{gather}
%
where $D(t)$ is slip at time $t$, $V$ is the slip rate, and $t_0$ is the slip initiation time (time when rupture reaches the location).
The slip rate is specified independently for each of the components of slip, and the slip rate and slip starting time may vary over the fault surface.

```{table} Values in the auxiliary field spatial database for KinSrcConstRate.
:name: tab:slip:function:step
|  Subfield         |  Components                          |
|:------------------|:-------------------------------------|
| `initiation_time` |            --                        |
| `slip_rate`       | `opening`, `left_lateral`, `reverse` |
```

