# Viewing fields

In addition to using the debugger to inspect code and variables, it is often helpful to print fields to stdout or inspect where a computed field does not match the expected field.
Turning on this type of output is usually done by activating debug journals via command line arguments.

Viewing a field will print the subfield metadata, the layout of the field, and the field values.
See {ref}`sec-developer-petsc-section` for how to interpret the layout of a field.

```{code-block} console
---
caption: Output from `Field::view()` for solution field with displacement and fluid pressure subfields.
---
Viewing field 'solution' Solution Field.
  Subfields: # Order of subfields is given by the index, not the order listed.
    Subfield displacement, index: 0, components: displacement_x displacement_y, scale: 1000, basisOrder: 1, quadOrder: 1
    Subfield fluid_pressure, index: 2, components: fluid_pressure, scale: 0.1, basisOrder: 1, quadOrder: 1
    Subfield velocity, index: 1, components: velocity_x velocity_y, scale: 100, basisOrder: 1, quadOrder: 1
  dimensionalize flag: 0
DM Object: 1 MPI processes
  type: plex
DM_0xe6a550_38 in 2 dimensions:
  0-cells: 5
  1-cells: 8
  2-cells: 4
Labels:
  boundary_bottom: 1 strata with value/size (1 (3))
  boundary: 1 strata with value/size (1 (8))
  material-id: 1 strata with value/size (24 (4))
  depth: 3 strata with value/size (0 (5), 1 (8), 2 (4))
PetscSection Object: 1 MPI processes
  type not yet set
3 fields
  field 0 with 2 components # displacement vector field
Process 0:
  (   0) dim  0 offset   0
  (   1) dim  0 offset   0
  (   2) dim  0 offset   0
  (   3) dim  0 offset   0
  (   4) dim  2 offset   0 constrained 1 # y degree of freedom is constrained
  (   5) dim  2 offset   5 constrained 1 # y degree of freedom is constrained
  (   6) dim  2 offset  10
  (   7) dim  2 offset  15
  (   8) dim  2 offset  20
  (   9) dim  0 offset  25
  (  10) dim  0 offset  25
  (  11) dim  0 offset  25
  (  12) dim  0 offset  25
  (  13) dim  0 offset  25
  (  14) dim  0 offset  25
  (  15) dim  0 offset  25
  (  16) dim  0 offset  25
  field 1 with 2 components # velocity vector field
Process 0:
  (   0) dim  0 offset   0
  (   1) dim  0 offset   0
  (   2) dim  0 offset   0
  (   3) dim  0 offset   0
  (   4) dim  2 offset   2
  (   5) dim  2 offset   7
  (   6) dim  2 offset  12
  (   7) dim  2 offset  17
  (   8) dim  2 offset  22
  (   9) dim  0 offset  25
  (  10) dim  0 offset  25
  (  11) dim  0 offset  25
  (  12) dim  0 offset  25
  (  13) dim  0 offset  25
  (  14) dim  0 offset  25
  (  15) dim  0 offset  25
  (  16) dim  0 offset  25
  field 2 with 1 components # pressure scalar field
Process 0:
  (   0) dim  0 offset   0
  (   1) dim  0 offset   0
  (   2) dim  0 offset   0
  (   3) dim  0 offset   0
  (   4) dim  1 offset   4
  (   5) dim  1 offset   9
  (   6) dim  1 offset  14
  (   7) dim  1 offset  19
  (   8) dim  1 offset  24
  (   9) dim  0 offset  25
  (  10) dim  0 offset  25
  (  11) dim  0 offset  25
  (  12) dim  0 offset  25
  (  13) dim  0 offset  25
  (  14) dim  0 offset  25
  (  15) dim  0 offset  25
  (  16) dim  0 offset  25
Proc 0 local vector # 5 nondimensionalized values per point: displacement (2), velocity (2), pressure (1)
Vec Object: unknown 1 MPI processes
  type: seq
-0.999 # offset  0, point 4, x-displacement
-4.2 # offset 1, point 4, y-displacement
-9.99 # offset 2, point 4, x-velocity
-9.99 # offset 3, point 4, y-velocity
-9990. # offset 4, point 4, pressure
-0.999 # offset 5, point 5, x-displacement
0.6
-9.99
-9.99
-9990.
-0.999 # offset 10, point 6, x-displacement
0.
-9.99
-9.99
-9990.
-0.999 # offset 15, point 7, x-displacement
-0.6
-9.99
-9.99
-9990.
-0.999 # offset 20, point 8, x-displacement
4.2
-9.99
-9.99
-9990.
```

## Viewing differences

In tests in which we compare a computed field against one from an analytical solution using
`DMPlexComputeL2DiffLocal()` and the fields do not agree, it is generally helpful to determine which pieces do not agree.
The `DMPlex` object contains an internal switch to print the point-by-point differences while computing the norm.
This switch can be activated using the `--petsc dm_plex_print_l2=1` command line argument in the C++ and MMS tests.

```{code-block} console
---
caption: Debugging output from `DMPlexComputeL2DiffLocal()`
---
Cell 0 Element Solution for Field 0 # displacement vector field
  | -0.999 | # Values of solution field variable at vertices of cell 0
  | -4.2 |
  | -0.999 |
  | 0. |
  | -0.999 |
  | -0.6 |
    elem 0 field 0 diff 0.  # Differences at quadrature points with respect to field given by analytical function
    elem 0 field 0 diff 6.27226e-32
    elem 0 field 0 diff 2.24281e-33
    elem 0 field 0 diff 8.97125e-33
    elem 0 field 0 diff 0.
    elem 0 field 0 diff 0.
    elem 0 field 0 diff 2.24281e-33
    elem 0 field 0 diff 0.
Cell 0 Element Solution for Field 1 # velocity field vector field
  | -9.99 |
  | -9.99 |
  | -9.99 |
  | -9.99 |
  | -9.99 |
  | -9.99 |
    elem 0 field 1 diff 0.
    elem 0 field 1 diff 0.
    elem 0 field 1 diff 0.
    elem 0 field 1 diff 0.
    elem 0 field 1 diff 0.
    elem 0 field 1 diff 0.
    elem 0 field 1 diff 0.
    elem 0 field 1 diff 0.
Cell 0 Element Solution for Field 2 # pressure scalar field
  | -9990. |
  | -9990. |
  | -9990. |
    elem 0 field 2 diff 0.
    elem 0 field 2 diff 0.
    elem 0 field 2 diff 0.
    elem 0 field 2 diff 0.
  elem 0 diff 7.61795e-32
Cell 1 Element Solution for Field 0
  | -0.999 | # Values of solution field variable at vertices of cell 1
  | 0.6 |
  | -0.999 |
  | 0. |
  | -0.999 |
  | -4.2 |
    elem 1 field 0 diff 0.
    elem 1 field 0 diff 3.92016e-33
    elem 1 field 0 diff 2.24281e-33
    elem 1 field 0 diff 1.4354e-31
    elem 1 field 0 diff 0.
    elem 1 field 0 diff 0.
    elem 1 field 0 diff 2.24281e-33
    elem 1 field 0 diff 0.
Cell 1 Element Solution for Field 1
  | -9.99 |
  | -9.99 |
  | -9.99 |
  | -9.99 |
  | -9.99 |
  | -9.99 |
    elem 1 field 1 diff 0.
    elem 1 field 1 diff 0.
    elem 1 field 1 diff 0.
    elem 1 field 1 diff 0.
    elem 1 field 1 diff 0.
    elem 1 field 1 diff 0.
    elem 1 field 1 diff 0.
    elem 1 field 1 diff 0.
Cell 1 Element Solution for Field 2
  | -9990. |
  | -9990. |
  | -9990. |
    elem 1 field 2 diff 0.
    elem 1 field 2 diff 0.
    elem 1 field 2 diff 0.
    elem 1 field 2 diff 0.
  elem 1 diff 1.51946e-31
... # Output continues for values in other cells
```
