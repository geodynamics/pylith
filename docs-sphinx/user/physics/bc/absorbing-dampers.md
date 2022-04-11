(sec-user-physics-absorbing-dampers)=
# Absorbing Boundary Conditions

You can use the `AbsorbingDampers` boundary condition to prevent seismic waves reflecting off of a boundary.
Normally incident dilatational and shear waves are perfectly absorbed.
Waves incident at other angles are only partially absorbed.
This boundary condition is simpler than a perfectly matched layer (PML) boundary condition but does not perform quite as well, especially for surface waves.
If the waves arriving at the absorbing boundary are relatively small in amplitude compared to the amplitudes of primary interest, this boundary condition gives reasonable results.

The auxiliary field spatial database contains the bulk rheology properties for an isotrpoic, linear elastic material (density, Vs (S-wave speed), and Vp (P-wave speed).
You can simply use the same spatial database that was used to specify the elastic properties of the material.

:::{seealso}
See [`AbsorbingDampers` Component](../../components/bc/AbsorbingDampers.md) for the Pyre properties and facilities and configuration examples.

See {ref}`sec-user-governing-eqns-absorbing-boundary` for the derivation of the finite-element implementation of the absorbing boundary.
:::
