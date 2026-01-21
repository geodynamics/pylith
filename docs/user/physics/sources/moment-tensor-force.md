(sec-user-physics-moment-tensor-force)=
# Moment Tensor Force (`MomentTensorForce`)

The `MomentTensorForce` component implements point sources using a moment tensor representation.
This is commonly used in seismology to represent earthquake sources or other internal force sources.

A moment tensor describes the equivalent body forces of a seismic source.
In 2D, the moment tensor has 4 independent components (symmetric $2 \times 2$ tensor), while in 3D it has 6 independent components (symmetric $3 \times 3$ tensor).

The temporal evolution of the source is controlled by a source time function, which multiplies the moment tensor amplitude at each time step.
PyLith provides several built-in source time functions as well as support for user-defined time histories.

## Source Parameters

The source parameters include the moment tensor components and a time delay.
The time delay allows offsetting the source time function for each point source, enabling simulation of rupture propagation.

```{table} Values in the auxiliary field spatial database for MomentTensorForce.
:name: tab:source:momenttensorforce:auxiliary
| Subfield           | Components (2D)            | Components (3D)                      |
|:-------------------|:---------------------------|:-------------------------------------|
| `moment_tensor`    | xx, yy, xy, zz             | xx, yy, zz, xy, yz, xz               |
| `time_delay`       | --                         | --                                   |
```

## Source Points File

The locations of point sources are specified in a text file.
The file contains the coordinates of each source point along with an optional name identifier.

```{code-block} text
---
caption: Example source points file format.
---
# Comment lines start with #
# x  y  [z]  name
0.0  0.0  source1
1000.0  500.0  source2
```

:::{seealso}
See [`MomentTensorForce` Component](../../components/sources/MomentTensorForce.md) for the Pyre properties and facilities and configuration examples.
:::

## Source Time Functions

The source time function controls the temporal evolution of the moment tensor force.
PyLith supports the following source time functions:

:`RickerWavelet`: A Ricker wavelet (Mexican hat wavelet) commonly used in seismic modeling.
:`GaussianWavelet`: A Gaussian wavelet.
:`SquareWavelet`: A step function (Heaviside function).
:`TimeHistoryWavelet`: A user-specified time function from a time history database.

### Ricker Wavelet (`RickerWavelet`)

The Ricker wavelet, also known as the Mexican hat wavelet, is the negative normalized second derivative of a Gaussian function.
It is commonly used as a source time function in seismic wave propagation studies because it has a well-defined frequency content and compact time duration.

The Ricker wavelet is defined as:
%
\begin{gather}
  S(t) = \left(1 - 2\pi^2 f_0^2 (t-t_d)^2\right) \exp\left(-\pi^2 f_0^2 (t-t_d)^2\right)
\end{gather}
%
where $f_0$ is the center frequency and $t_d$ is the time delay.

```{table} Values in the auxiliary field spatial database for RickerWavelet.
:name: tab:source:rickerwavelet:auxiliary
| Subfield            | Description                              |
|:--------------------|:-----------------------------------------|
| `center_frequency`  | Center frequency of the wavelet ($f_0$)  |
```

:::{seealso}
See [`RickerWavelet` Component](../../components/sources/RickerWavelet.md) for the Pyre properties and facilities and configuration examples.
:::

### Gaussian Wavelet (`GaussianWavelet`)

The Gaussian wavelet provides a smooth pulse.
It is defined as:
%
\begin{gather}
  S(t) = \frac{1}{2\pi^2 f_0^2} \exp\left(-\pi^2 f_0^2 (t-t_d)^2\right)
\end{gather}
%
where $f_0$ is the center frequency and $t_d$ is the time delay.

:::{note}
The implementation uses the negative of the second derivative of a Gaussian.
:::

```{table} Values in the auxiliary field spatial database for GaussianWavelet.
:name: tab:source:gaussianwavelet:auxiliary
| Subfield            | Description                              |
|:--------------------|:-----------------------------------------|
| `center_frequency`  | Center frequency of the wavelet ($f_0$)  |
```

:::{seealso}
See [`GaussianWavelet` Component](../../components/sources/GaussianWavelet.md) for the Pyre properties and facilities and configuration examples.
:::

### Square Wavelet (`SquareWavelet`)

The square wavelet is a step function (Heaviside function) that activates the source at the time delay:
%
\begin{gather}
  S(t) = \left\{ \begin{array}{cc}
    0 & t < t_d \\
    1 & t \geq t_d
  \end{array}\right.
\end{gather}
%
where $t_d$ is the time delay.

```{table} Values in the auxiliary field spatial database for SquareWavelet.
:name: tab:source:squarewavelet:auxiliary
| Subfield            | Description                              |
|:--------------------|:-----------------------------------------|
| `center_frequency`  | Not used (placeholder for consistency)   |
```

:::{seealso}
See [`SquareWavelet` Component](../../components/sources/SquareWavelet.md) for the Pyre properties and facilities and configuration examples.
:::

### Time History Wavelet (`TimeHistoryWavelet`)

The time history wavelet allows specification of an arbitrary source time function using a time history database.
The time history should contain normalized amplitudes (typically between 0 and 1) as a function of time.

```{table} Values in the auxiliary field spatial database for TimeHistoryWavelet.
:name: tab:source:timehistorywavelet:auxiliary
| Subfield                   | Description                                   |
|:---------------------------|:----------------------------------------------|
| `time_history_start_time`  | Start time for the time history               |
```

:::{seealso}
See [`TimeHistoryWavelet` Component](../../components/sources/TimeHistoryWavelet.md) for the Pyre properties and facilities and configuration examples.
:::

## Example Configuration

The following example shows how to configure a moment tensor point source with a Ricker wavelet:

```{code-block} cfg
---
caption: Example configuration for a moment tensor point source with Ricker wavelet.
---
[pylithapp.problem]
sources = [source]
sources.source = pylith.sources.MomentTensorForce

[pylithapp.problem.sources.source]
description = Earthquake source
label_value = 2

# Specify source locations
reader.filename = source_sites.txt
reader.coordsys = spatialdata.geocoords.CSCart
reader.coordsys.space_dim = 2

# Use Ricker wavelet source time function
source_time_function = pylith.sources.RickerWavelet

# Source parameters (moment tensor and time delay)
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Source properties
db_auxiliary_field.values = [moment_tensor_xx, moment_tensor_yy, moment_tensor_xy, moment_tensor_zz, time_delay, center_frequency]
db_auxiliary_field.data = [1.0e12*Pa*s, 1.0e12*Pa*s, 0.0*Pa*s, 1.0e12*Pa*s, 0.0*s, 5.0]

# Discretization of auxiliary subfields
auxiliary_subfields.moment_tensor.basis_order = 0
auxiliary_subfields.time_delay.basis_order = 0
source_time_function.auxiliary_subfields.center_frequency.basis_order = 0
```

