# OutputTriggerTime

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.OutputTriggerTime`
:Journal name: `outputtriggertime`

Define how often output is written in terms of elasped simulation time.

:::{tip}
Due to floating point roundoff, it is usually a good idea to use a value that is a fraction of a time step smaller than the desired value.
:::

Implements `OutputTrigger`.

## Pyre Properties

* `elapsed_time`=\<dimensional\>: Elapsed time between writes.
  - **default value**: 0*s
  - **current value**: 0*s, from {default}

## Example

Example of setting `OutputTriggerTime` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[output_trigger]
elapsed_time = 0.9999*year
:::

