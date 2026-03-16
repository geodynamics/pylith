# MomentTensorForce

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.sources.MomentTensorForce`
:Journal name: `momenttensorforce`

Moment tensor point source.

This component implements a point source using a moment tensor representation.
The moment tensor describes the equivalent body forces of a seismic source.
The temporal evolution of the source is controlled by a source time function.

Implements `Source`.

## Pyre Facilities

* `auxiliary_subfields`: Discretization information for auxiliary subfields.
  - **current value**: 'auxiliary_subfields', from {file='...', line=..., function='__set__'}
  - **configurable as**: auxiliary_subfields
* `db_auxiliary_field`: Database for source parameters.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db_auxiliary_field
* `derived_subfields`: Discretization of derived subfields.
  - **current value**: 'emptybin', from {default}
  - **configurable as**: emptybin, derived_subfields
* `observers`: Observers (e.g., output).
  - **current value**: 'singlephysicsobserver', from {default}
  - **configurable as**: singlephysicsobserver, observers
* `reader`: Reader for source points list.
  - **current value**: 'pointslist', from {default}
  - **configurable as**: pointslist, reader
* `source_time_function`: Source time function for moment tensor force.
  - **current value**: 'timehistorywavelet', from {default}
  - **configurable as**: timehistorywavelet, source_time_function

## Pyre Properties

* `description`=\<str\>: Descriptive label for source.
  - **default value**: ''
  - **current value**: '', from {default}
  - **validator**: <function validateDescription at ...>
* `field`=\<str\>: Solution subfield associated with source.
  - **default value**: 'displacement'
  - **current value**: 'displacement', from {default}
* `label`=\<str\>: Name of label for source.
  - **default value**: 'source-id'
  - **current value**: 'source-id', from {default}
* `label_value`=\<int\>: Value of label identifying source.
  - **default value**: 1
  - **current value**: 1, from {default}

## Example

Example of setting `MomentTensorForce` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
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
:::

