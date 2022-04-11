(sec-user-physics-slip-impulses)=
# Green's Functions (`FaultCohesiveImpulses`)

Computing static Green's functions using the `GreensFns` problem ({ref}`sec-user-problems-greensfns`) requires a specialized fault implementation, `FaultCohesiveImpulses`, to set up the slip impulses.
The parameters controlling the slip impulses include the slip components to be used (left-lateral, reverse, or fault opening) and the amplitude of the impulses (for example, to use only a subset of a fault).

:::{admonition} TODO
:class: error

Add information about slip discretization (0=cell, 1=vertices, etc).
:::

:::{seealso}
See [`FaultCohesiveImpulses` Component](../../components/faults/FaultCohesiveImpulses.md)  for the Pyre properties and facilities and configuration examples.
:::