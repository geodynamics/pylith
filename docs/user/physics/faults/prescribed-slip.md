(sec-user-physics-prescribed-slip)=
# Prescribed Slip (`FaultCohesiveKin`)

Prescribed slip, often called kinematic earthquake ruptures, use the `FaultCohesiveKin` component to prescribe the slip as a function of time on the fault surface.
Slip may evolve simultaneously over the fault surface instantaneously in a single time step (as is usually done in quasistatic simulations) or propagate over the fault surface over hundreds and up to thousands of time steps (as is usually done in a dynamic simulation).

Multiple earthquake ruptures can be specified on a single fault surface.
This permits repeatedly rupturing the same portion of a fault or combining earthquake rupture on one subset of the fault surface with steady aseismic slip on another subset (the two subsets may overlap in both time and space).
An array of kinematic earthquake rupture components associates a name (string) with each kinematic rupture.
The default dynamic array contains a single earthquake rupture, `rupture`.

The default discretization of the `slip` auxiliary subfield is a basis order of 1 (linear variation on the fault).
This discretization will also be used for all auxiliary subfields associated with parameters of the slip time function.

:::{seealso}
See [`FaultCohesiveKin` Component](../components/faults/../../../components/faults/FaultCohesiveKin.md) for the Pyre properties and facilities and configuration examples.
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
    0 & 0\leq t < t_r  \\
    D_{final} & t\ge t_r
  \end{array}\right.
\end{gather}
%
where $D(t)$ is slip at time $t$, $D_{final}$ is the final slip, and $t_r$ is the slip initiation time (time when rupture reaches the location).
The slip is specified independently for each of the components of slip, and the slip and slip starting time may vary over the fault surface.

:::{figure-md} fig:slipfn:step
<img src="figs/slipfn-step.*" alt="Step slip time function." width="600px"/>

Step slip time function.
:::

```{table} Values in the auxiliary field spatial database for KinSrcStep.
:name: tab:slip:function:step
|      Subfield     |           Components                 |
|:------------------|:-------------------------------------|
| `initiation_time` ($t_r$) |     --                               |
| `final_slip`      | `opening`, `left_lateral`, `reverse` |
```

:::{seealso}
See [`KinSrcStep` Component](../components/faults/../../../components/faults/KinSrcStep.md) for the Pyre properties and facilities and configuration examples.
:::

### Constant Slip Rate Slip Time Function (`KinSrcConstRate`)

This slip function prescribes a constant slip rate for the evolution of slip at a point:
%
\begin{gather}
  D(t)=\left\{ \begin{array}{cc}
    0 & 0\leq t < t_r \\
    V(t-t_r) & t \ge t_r
  \end{array}\right.
\end{gather}
%
where $D(t)$ is slip at time $t$, $V$ is the slip rate, and $t_r$ is the slip initiation time (time when rupture reaches the location).
The slip rate is specified independently for each of the components of slip, and the slip rate and slip starting time may vary over the fault surface.

:::{figure-md} fig:slipfn:constrate
<img src="figs/slipfn-constrate.*" alt="Constant rate slip time function." width="600px"/>

Constant rate slip time function.
:::

```{table} Values in the auxiliary field spatial database for KinSrcConstRate.
:name: tab:slip:function:constrate
|  Subfield         |  Components                          |
|:------------------|:-------------------------------------|
| `initiation_time` ($t_r$) |            --                        |
| `slip_rate`       | `opening`, `left_lateral`, `reverse` |
```

:::{seealso}
See [`KinSrcConstRate` Component](../components/faults/../../../components/faults/KinSrcConstRate.md) for the Pyre properties and facilities and configuration examples.
:::

### Ramp Slip Time Function (`KinSrcRamp`)

This slip function prescribes a constant slip rate over a time window with a smooth initiation and termination:
%
\begin{gather}
  D(t)=\left\{ \begin{array}{cc}
%  
    0 & t-t_r \leq 0 \\
%
  \frac{D_\mathit{fina}}{C_\mathit{acc}} f_1(t) & t-t_r \leq \frac{1}{2} t_\mathit{acc} \\
  \frac{D_\mathit{fina}}{C_\mathit{acc}} f_2(t) & t-t_r \leq t_\mathit{acc} \\
  \frac{D_\mathit{fina}}{C_\mathit{acc}} f_3(t) & t-t_r \leq t_\mathit{rise}-t_\mathit{acc}) \\
  \frac{D_\mathit{fina}}{C_\mathit{acc}} f_4(t) & t-t_r \leq t_\mathit{rise}-\frac{1}{2} t_\mathit{acc} \\
  \frac{D_\mathit{fina}}{C_\mathit{acc}} f_5(t) & t-t_r \leq t_\mathit{rise} \\
  \frac{D_\mathit{fina}}{C_\mathit{acc}} f_6(t) & t-t_r > t_\mathit{rise}
\end{array}\right.
\end{gather}

\begin{align}
f_1(t) = &\frac{1}{6} t^3 \\
%
f_2(t) =
    &\frac{1}{2} t_\mathit{acc} t^2
    -\frac{1}{6} t^3
    -\frac{1}{4} t_\mathit{acc}^2 t
    +\frac{1}{24} t_\mathit{acc}^3 \\
%
f_3(t) =
    &\frac{1}{4} t_\mathit{acc}^2 t
    -\frac{1}{8} t_\mathit{acc}^3 \\
%
f_4(t) =
    &-\frac{1}{6} t^3
    +\frac{1}{2} (t_\mathit{rise}-t_\mathit{acc}) t^2
    -\frac{1}{2} (t_\mathit{rise}-t_\mathit{acc})^2 t
    +\frac{1}{4} t_\mathit{acc}^2 t
    +\frac{1}{6} (t_\mathit{rise}-t_\mathit{acc})^3
    -\frac{1}{8} t_\mathit{acc}^3 \\
%
f_5(t) =
    &\frac{1}{6} t^3
    -\frac{1}{2} t_\mathit{rise} t^2
    +\frac{1}{2} t_\mathit{rise}^2 t
    -\frac{1}{3} t_1^3
    +\frac{1}{2} t_\mathit{rise} t_1^2
    -\frac{1}{2} t_\mathit{rise}^2 t_1 \\
    &+\frac{1}{2} (t_\mathit{rise}-t_\mathit{acc}) t_1^2
    -\frac{1}{2} (t_\mathit{rise}-t_\mathit{acc})^2 t_1
    +\frac{1}{4} t_\mathit{acc}^2  t_1
    +\frac{1}{6} (t_\mathit{rise}-t_\mathit{acc})^3
    -\frac{1}{8} t_\mathit{acc}^3 \\
%
f_6(t) =
    &\frac{1}{6} (t_\mathit{rise}^3
    -\frac{1}{3} (t_1^3
    +\frac{1}{2} t_\mathit{rise} (t_1^2
    -\frac{1}{2} (t_\mathit{rise}^2 t_1
    +\frac{1}{2} (t_\mathit{rise}-t_\mathit{acc}) t_1^2
    -\frac{1}{2} (t_\mathit{rise}-t_\mathit{acc})^2 t_1 \\
    &+\frac{1}{4} t_\mathit{acc}^2 t_1
    +\frac{1}{6} (t_\mathit{rise}-t_\mathit{acc})^3
    -\frac{1}{8} t_\mathit{acc}^3 \\
C_\mathit{acc} = 
    &\frac{1}{6} t_\mathit{rise}^3
    -1.0/3.0 (t_\mathit{rise}-\frac{1}{2} t_\mathit{acc})^3
    +\frac{1}{2} t_\mathit{rise} (t_\mathit{rise}-\frac{1}{2} t_\mathit{acc})^2
    -\frac{1}{2} t_\mathit{rise}^2  (t_\mathit{rise}-\frac{1}{2} t_\mathit{acc}) \\
    &+\frac{1}{2} (t_\mathit{rise}-t_\mathit{acc}) (t_\mathit{rise}-\frac{1}{2} t_\mathit{acc})^2
    -\frac{1}{2} (t_\mathit{rise}-t_\mathit{acc})^2 (t_\mathit{rise}-\frac{1}{2} t_\mathit{acc})
    +\frac{1}{4} (t_\mathit{acc}^2 (t_\mathit{rise}-\frac{1}{2} t_\mathit{acc})\\
    &+\frac{1}{6} (t_\mathit{rise}-t_\mathit{acc})^3
    -\frac{1}{8} (t_\mathit{acc})^3 \\
%
t_1 = &t_\mathit{rise} - \frac{1}{2} t_\mathit{acc}
\end{align}
%
where $D(t)$ is slip at time $t$, $t_r$ is the slip initiation time (time when rupture reaches the location), $t_\mathit{rise}$ is the rise time, and $t_\mathit{acc}$ is the duration of the acceleration impulse.
The slip is specified independently for each of the components of slip, and the slip rate and slip starting time may vary over the fault surface.

:::{figure-md} fig:slipfn:ramp
<img src="figs/slipfn-ramp.*" alt="Ramp slip time function." width="600px"/>

Ramp slip time function.
:::

```{table} Values in the auxiliary field spatial database for KinSrcRamp.
:name: tab:slip:function:ramp
|  Subfield          |  Components                          |
|:-------------------|:-------------------------------------|
| `initiation_time` ($t_r$)  |            --                        |
| `rise_time` ($t_\mathit{rise}$)       |            --                        |
| `impulse_duration` ($t_\mathit{acc}$) |            --                        |
| `final_slip`       | `opening`, `left_lateral`, `reverse` |
```

:::{seealso}
See [`KinSrcRamp` Component](../components/faults/../../../components/faults/KinSrcRamp.md) for the Pyre properties and facilities and configuration examples.
:::

### Brune Slip Time Function (`KinSrcBrune`)

This slip function corresponds to the integral of Brune's far-field time function {cite}`Brune:1970` and prescribes a rapid rise in slip rate followed by a very gradual slowdown:
%
\begin{gather}
  D(t) = \left\{ \begin{array}{cc}
    0 & 0\leq t < t_r \\
    D_\mathit{final} \left(1-exp\left(-\frac{t-t_r}{t_1}\right)\left(1+\frac{t-t_r}{t_1}\right)\right) & t \ge t_r
  \end{array}\right.\\
  t_1 = 0.6195 t_\mathit{rise}
\end{gather}
where $D(t)$ is slip at time $t$, $D_{final}$ is the final slip at the location, $t_r$ is the slip initiation time (time when rupture reaches the location), and $t_\mathit{rise}$ is the rise time.
Because the slip time function approaches the final slip asymptotically, we use the time it takes for the slip to reach 95\% of the final slip value as the rise time.

:::{figure-md} fig:slipfn:brune
<img src="figs/slipfn-brune.*" alt="Brune slip time function." width="600px"/>

Brune slip time function.
:::

```{table} Values in the auxiliary field spatial database for KinSrcBrune.
:name: tab:slip:function:brune
|  Subfield                       |  Components                          |
|:--------------------------------|:-------------------------------------|
| `initiation_time` ($t_r$)       |            --                        |
| `rise_time` ($t_\mathit{rise}$) |            --                        |
| `final_slip`                    | `opening`, `left_lateral`, `reverse` |
```

:::{seealso}
See [`KinSrcBrune` Component](../components/faults/../../../components/faults/KinSrcBrune.md) for the Pyre properties and facilities and configuration examples.
:::

### Liu-Cosine Slip Time Function (`KinSrcLiuCosine`)

This slip time function, proposed by Liu, Archuleta, and Hartzell for use in ground-motion modeling {cite}`Liu:etal:2006`, combines several cosine and sine functions together to create a slip time history with a sharp rise and gradual termination with a finite duration of slip.
The evolution of slip at a point follows:
%
\begin{gather}
  D(t) = \left\{ \begin{array}{cc}
    D_{\mathit{final}} C_n \left(0.7t-0.7\frac{t_1}\pi\sin\frac{\pi t}{t_1}-1.2\frac{t_1}{\pi}\left(\cos\frac{\pi t}{2t_1}-1\right)\right) & 0\leq t < t_1 \\
%
    D_{\mathit{final}} C_n \left(1.0 t - 0.7\frac{t1}{\pi}\sin\frac{\pi t}{t_1} + 0.3\frac{t2}{\pi}\sin\frac{\pi(t-t1)}{t_2}+\frac{1.2}{\pi}t_1-0.3t_1\right) & t_1\leq t < 2t_1\\
%
    D_{\mathit{final}}C_n\left(0.7-0.7\cos\frac{\pi t}{t_1}+0.6\sin\frac{\pi t}{2 t_1}\right) & 2 t_1 \leq t \leq t_0
  \end{array}\right.\\
  C_n=\frac{\pi}{1.4\pi t_1+1.2t_1+0.3\pi t_2}\\
  t_0=1.525 t_\mathit{rise}\\
  t_1=0.13t_0\\
  t_2=t_0-t_1
\end{gather}
%
where $D(t)$ is slip at time $t$, $D_{final}$ is the final slip at the location, $t_r$ is the slip initiation time (time when rupture reaches the location), and $t_\mathit{rise}$ is the rise time.

:::{figure-md} fig:slipfn:liucos
<img src="figs/slipfn-liucos.*" alt="Liu-cosine slip time function." width="600px"/>

Liu-cosine slip time function.
:::

```{table} Values in the auxiliary field spatial database for KinSrcLiuCos.
:name: tab:slip:function:liucos
|  Subfield                       |  Components                          |
|:--------------------------------|:-------------------------------------|
| `initiation_time` ($t_r$)       |            --                        |
| `rise_time` ($t_\mathit{rise}$) |            --                        |
| `final_slip`                    | `opening`, `left_lateral`, `reverse` |
```

:::{seealso}
See [`KinSrcLiuCos` Component](../components/faults/../../../components/faults/KinSrcLiuCos.md) for the Pyre properties and facilities and configuration examples.
:::

### User-Time History Slip Time Function (`KinSrcTimeHistory`)

This slip time function reads the slip time function from a data file, so it can have an arbitrary shape.
The slip and slip initiation times are specified using spatial databases, so the slip time function should use a normalized amplitude (0 $\rightarrow$ 1).

```{table} Values in the auxiliary field spatial database for KinSrcTimeHistory.
:name: tab:slip:function:timehistory
|  Subfield         |  Components                          |
|:------------------|:-------------------------------------|
| `initiation_time` |            --                        |
| `final_slip`      | `opening`, `left_lateral`, `reverse` |
```

:::{seealso}
See [`KinSrcTimeHistory` Component](../components/faults/../../../components/faults/KinSrcTimeHistory.md) for the Pyre properties and facilities and configuration examples.
:::
