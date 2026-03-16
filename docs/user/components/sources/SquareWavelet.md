# SquareWavelet

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.sources.SquareWavelet`
:Journal name: `squarewavelet`

Square wavelet (step function) source time function.

The square wavelet is a Heaviside step function that activates the source at the time delay.

The wavelet is defined as:

$S(t) = \begin{cases} 0 & t < t_d \\ 1 & t \geq t_d \end{cases}$

where $t_d$ is the time delay.

Implements `SourceTimeFunctionMomentTensorForce`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for source time function parameters.
  - **current value**: 'auxiliary_subfields', from {file='...', line=..., function='__set__'}
  - **configurable as**: auxiliary_subfields

## Example

Example of setting `SquareWavelet` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp.problem.sources.source]
source_time_function = pylith.sources.SquareWavelet

# Source parameters specified in the auxiliary field database
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Source properties
db_auxiliary_field.values = [moment_tensor_xx, moment_tensor_yy, moment_tensor_xy, moment_tensor_zz, time_delay, center_frequency]
db_auxiliary_field.data = [1.0e12*Pa*s, 1.0e12*Pa*s, 0.0*Pa*s, 1.0e12*Pa*s, 0.5*s, 1.0]

source_time_function.auxiliary_subfields.center_frequency.basis_order = 0
:::

