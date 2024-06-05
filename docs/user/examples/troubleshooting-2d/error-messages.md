# Error Messages

Summary of error messages addressed in this set of examples.

## Configuration Errors

```{code-block} bash
---
caption: See [Step 1, Error 1](step01-error01.md).
---
 >> {default}::
 -- pyre.inventory(error)
 -- timedependent.problem_defaults.name <- ''
 -- Missing required property 'name' in default options for problem.
```

```{code-block} bash
---
caption: See [Step 1, Error 2](step01-error02.md).
---
 >> {default}::
 -- pyre.inventory(error)
 -- timedependent.bc.dirichlettimedependent.label <- ''
 -- Label for boundary condition group/nodeset/pset in mesh not specified.
```

```{code-block} bash
---
caption: See [Step 1, Error 3](step01-error03.md).
---
>> ./pylithapp.cfg:133:
-- pyre.inventory(error)
-- pylithapp.timedependent.materials.elasticity.auxiliary_subfields.bulk_modulus.basis_order <- '0'
-- unknown component 'pylithapp.timedependent.materials.elasticity.auxiliary_subfields.bulk_modulus'
```

```{code-block} bash
---
caption: See [Step 1, Error 4](step01-error04.md).
---
RuntimeError: Material label_value '3' for cell '609' does not match the label_value of any materials or interfaces.
```

```{code-block}
---
caption: See [Step 6, Error 1](step06-error01.md).
---
>> {default}::
-- pyre.inventory(error)
-- timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb.description <- ''
-- Description for spatial database not specified.
```

```{code-block}
---
caption: See [Step 6, Error 2](step06-error02.md).
---
>> {default}::
-- pyre.inventory(error)
-- timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb.simpleioascii.filename <- ''
-- Filename for spatial database not specified.
```

```{code-block}
---
caption: See [Step 6, Error 3](step06-error03.md).
---
RuntimeError: Cannot find 'lagrange_multiplier_fault' subfield in solution field for fault implementation in component 'splay'.
```

## Runtime Errors

```{code-block} bash
---
caption: See [Step 1, Error 5](step01-error05.md).
---
[0]PETSC ERROR: Error converting spatial database values for gravitational_acceleration at (  -96623.5  -72650.4) in
spatial database 'Gravity field'. Found near zero magnitude (0) for gravity field vector (  0  0).
```

```{code-block} bash
---
caption: See [Step 6, Error 4](step06-error04.md).
---
RuntimeError: Error occurred while reading spatial database file 'fault_slip.spatialdb'.
Read data for 3 out of 4 points.
```

```{code-block} bash
---
caption: See [Step 6, Error 5](step06-error05.md).
---
RuntimeError: Could not find value 'final_slip_opening' in spatial database 'Fault rupture for main fault'. Available values are:
  final-slip-left-lateral
  final-slip-opening
  initiation-time
```

```{code-block} bash
---
caption: See [Step 6, Error 6](step06-error06.md).
---
Nonlinear solve did not converge due to DIVERGED_MAX_IT iterations 1
```

```{code-block} bash
---
caption: See [Step 6, Error 7](step06-error07.md).
---
[0]PETSC ERROR: Residual norm computed by GMRES recursion formula 3.48613e+10 is far from the computed residual norm 6.92443e+12 at restart, residual norm at start of cycle 6.91369e+12
```

```{code-block} bash
---
caption: See [Step 6, Error 9](step06-error09.md).
---
[0]PETSC ERROR: Could not find values for initiation_time at (  -21860.3  -27621) in spatial database 'Fault rupture for main fault'.
```
