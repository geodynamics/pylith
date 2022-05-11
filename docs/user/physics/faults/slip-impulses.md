(sec-user-physics-fault-cohesive-impulses)=
# Fault Slip Impulses (`FaultCohesiveImpulses`)

Computing static Green's functions using the `GreensFns` problem ({ref}`sec-user-problems-greensfns`) requires a specialized fault implementation, `FaultCohesiveImpulses`, to set up the slip impulses.
The parameters controlling the slip impulses include the slip components to be used (left-lateral, reverse, or fault opening) and the amplitude of the impulses (for example, to use only a subset of a fault).

The slip auxiliary subfield specifies the patch of the fault on which to apply slip impulses.
Additionally, the basis order for the representation of the slip auxiliary subfield and the solution fields controls the spatial representation of the slip impulse.
The default basis order for all of these fields is 1, corresponding to a linear variation in slip, fault traction, and displacement.
A basis order of 0 will approximate uniform slip within a fault cell using the finite-element representation of the displacement field and fault Lagrange multiplier.
For example, a basis order of 0 for the slip auxiliary subfield and a basis order of 1 for the displacement and fault Lagrange multiplier fields will result in nearly uniform slip within the fault cell with a linear taper into the adjacent surrounding fault cells.

Impulses will be applied at any point on the fault with a slip component greater than the specified threshold. Slip impulses can be applied on any subset of the slip components, such as just the left lateral and reverse components in 3D and not fault opening.

```{table} Impulse degrees of freedom
:name: tab:fault:impulses:impulse:dof
| Impulse degree of freedom | Slip component    |
| :------------------------ | :---------------- |
| 0                         | fault opening     |
| 1                         | left lateral      |
| 2                         | reverse (3D only) |
```

:::{admonition} Pyre User Interface
:class: seealso
See [`FaultCohesiveImpulses` Component](../../components/faults/FaultCohesiveImpulses.md)  for the Pyre properties and facilities and configuration examples.
:::
