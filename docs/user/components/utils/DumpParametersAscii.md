# DumpParametersAscii

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.utils.DumpParametersAscii`
:Journal name: `dumpparamters`

Dump PyLith parameter information to an ASCII file.

If you do not set the filename for the progress monitor, then PyLith will create one using the
simulation name from the application defaults settings.

Implements `DumpParameters`.

## Pyre Properties

* `filename`=\<str\>: Name of file written with parameters.
  - **default value**: ''
  - **current value**: '', from {default}
* `indent`=\<int\>: Nmber of spaces to indent.
  - **default value**: 4
  - **current value**: 4, from {default}
* `verbose`=\<bool\>: Include description, location, and aliases.
  - **default value**: True
  - **current value**: True, from {default}

## Example

Example of setting `DumpParametersAscii` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp]
dump_parameters = pylith.utils.DumpParametersAscii

[pylithapp.dump_parameters]
filename = output/parameters.txt
verbose = True
:::

