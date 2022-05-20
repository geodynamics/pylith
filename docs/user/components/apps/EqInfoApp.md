# EqInfoApp

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.apps.EqInfoApp`
:Journal name: `eqinfoapp`


## Pyre Facilities

* `coordsys`: Coordinate system associated with mesh.
  - **current value**: 'cscart', from {default}
  - **configurable as**: cscart, coordsys
* `db_properties`: Spatial database for elastic properties.
  - **current value**: 'simpledb', from {default}
  - **configurable as**: simpledb, db_properties
* `weaver`: the pretty printer of my configuration as an XML document
  - **current value**: 'weaver', from {default}
  - **configurable as**: weaver

## Pyre Properties

* `faults`=\<list\>: Array of fault names.
  - **default value**: []
  - **current value**: [], from {default}
* `filename_pattern`=\<str\>: Pattern for fault files.
  - **default value**: 'output/fault_%s.h5'
  - **current value**: 'output/fault_%s.h5', from {default}
* `output_filename`=\<str\>: Filename for output.
  - **default value**: 'eqstats.py'
  - **current value**: 'eqstats.py', from {default}
* `snapshot_units`=\<dimensional\>: Units for timestamps in array of snapshots.
  - **default value**: 1*s
  - **current value**: 1*s, from {default}
* `snapshots`=\<list\>: Array of timestamps for slip snapshots (-1 == last time step).
  - **default value**: [-1]
  - **current value**: [-1], from {default}
* `typos`=\<str\>: Specifies the handling of unknown properties and facilities
  - **default value**: 'pedantic'
  - **current value**: 'pedantic', from {default}
  - **validator**: (in ['relaxed', 'strict', 'pedantic'])

