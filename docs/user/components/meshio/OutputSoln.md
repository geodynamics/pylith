# OutputSoln

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.OutputSoln`
:Journal name: `outputsoln`

Abstract base class for output of solution subfields.

Implements `OutputObserver`.

## Pyre Facilities

* `trigger`: Trigger defining how often output is written.
  - **current value**: 'outputtriggerstep', from {default}
  - **configurable as**: outputtriggerstep, trigger
* `writer`: Writer for data.
  - **current value**: 'datawriterhdf5', from {default}
  - **configurable as**: datawriterhdf5, writer

## Pyre Properties

* `data_fields`=\<list\>: Names of solution subfields to include in output.
  - **default value**: ['all']
  - **current value**: ['all'], from {default}
* `output_basis_order`=\<int\>: Basis order for output.
  - **default value**: 1
  - **current value**: 1, from {default}
  - **validator**: (in [0, 1])

