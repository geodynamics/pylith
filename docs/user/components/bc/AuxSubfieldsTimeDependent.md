# AuxSubfieldsTimeDependent

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.bc.AuxSubfieldsTimeDependent`
:Journal name: `auxfieldstimedependent`

Auxiliary subfields for time-dependent boundary conditions.

The boundary conditions values have the functional form:

\begin{equation}
  f(x,t) = f_0(x) + \dot{f}_1(x)(t-t_1(x)) + f_2(x)a(t-t_2(x))
\end{equation}

The association of these functions with the auxiliary subfields is:

:$f_0(x)$: `initial_amplitude`
:$\dot{f}_1(x)$: `rate_amplitude`
:$t_1(x)$: `rate_start`
:$f_2(x)$: `time_history_amplitude`
:$t_2(x)$: `time_history_start`

## Pyre Facilities

* `initial_amplitude`: Initial amplitude, f_0(x), subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, initial_amplitude
* `rate_amplitude`: Rate amplitude, \dot{f}_1(x), subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, rate_amplitude
* `rate_start_time`: Rate starting time, t_1(x), subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, rate_start_time
* `time_history_amplitude`: Time history amplitude, f_2(x). subfield
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, time_history_amplitude
* `time_history_start_time`: Time history starting time, t_2(s), subfield.
  - **current value**: 'subfield', from {default}
  - **configurable as**: subfield, time_history_start_time

## Example

Example of setting `AuxSubfieldsTimeDependent` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[time_dependent_subfields]
initial_amplitude.basis_order = 1
rate_amplitude.basis_order = 0
rate_start_time.basis_order = 1
time_history_amplitude.basis_order = 1
time_history_start_time.basis_order = 1
:::

