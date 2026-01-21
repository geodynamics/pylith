# RickerWavelet

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.sources.RickerWavelet`
:Journal name: `rickerwavelet`

Ricker wavelet source time function.

The Ricker wavelet (Mexican hat wavelet) is the negative normalized second derivative of a Gaussian.
It is commonly used as a source time function in seismic wave propagation studies.

The wavelet is defined as:

$S(t) = \left(1 - 2\pi^2 f_0^2 (t-t_d)^2\right) \exp\left(-\pi^2 f_0^2 (t-t_d)^2\right)$

where $f_0$ is the center frequency and $t_d$ is the time delay.

Implements `SourceTimeFunctionMomentTensorForce`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for source time function parameters.
  - **current value**: 'auxiliary_subfields', from {file='...', line=..., function='__set__'}
  - **configurable as**: auxiliary_subfields

## Example

Example of setting `RickerWavelet` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.sources.source]
source_time_function = pylith.sources.RickerWavelet

# Center frequency specified in the auxiliary field database
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Source properties
db_auxiliary_field.values = [moment_tensor_xx, moment_tensor_yy, moment_tensor_xy, moment_tensor_zz, time_delay, center_frequency]
db_auxiliary_field.data = [1.0e12*Pa*s, 1.0e12*Pa*s, 0.0*Pa*s, 1.0e12*Pa*s, 0.0*s, 5.0]

source_time_function.auxiliary_subfields.center_frequency.basis_order = 0
:::

