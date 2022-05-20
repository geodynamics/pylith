# OutputTriggerStep

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.meshio.OutputTriggerStep`
:Journal name: `outputtriggerstep`

Define how often output is written in terms of solution steps.

Implements `OutputTrigger`.

## Pyre Properties

* `num_skip`=\<int\>: Number of solution steps to skip between writes (0 means write every time step).
  - **default value**: 0
  - **current value**: 0, from {default}
  - **validator**: (greater than or equal to 0)

## Example

Example of setting `OutputTriggerStep` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[output_trigger]
num_skip = 2
:::

