# Error Messages

Summary of error messages addressed in this set of examples.

## Configuration Errors

```{code-block} pyrejournal
---
caption: Refer to [Step 1, Error 1](step01-error01.md).
---
 >> {default}::
 -- error (pyre.inventory)
 -- timedependent.problem_defaults.name <- ''
 -- Missing required property 'name' in default options for problem.
```

```{code-block} pyrejournal
---
caption: Refer to [Step 1, Error 2](step01-error02.md).
---
 >> {default}::
 -- error (pyre.inventory)
 -- timedependent.bc.dirichlettimedependent.label <- ''
 -- Label for boundary condition group/nodeset/pset in mesh not specified.
```

```{code-block} pyrejournal
---
caption: Refer to [Step 1, Error 3](step01-error03.md).
---
 >> ./pylithapp.cfg:112:
 -- error (pyre.inventory)
 -- pylithapp.timedependent.materials.elasticity.auxiliary_subfields.bulk_modulus.basis_order <- '0'
 -- unknown component 'pylithapp.timedependent.materials.elasticity.auxiliary_subfields.bulk_modulus'
```

```{code-block} pyrejournal
---
caption: Refer to [Step 6, Error 1](step06-error01.md).
---
 >> {default}::
 -- error (pyre.inventory)
 -- timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb.description <- ''
 -- Description for spatial database not specified.
```

```{code-block} pyrejournal
---
caption: Refer to [Step 6, Error 2](step06-error02.md).
---
 >> {default}::
 -- error (pyre.inventory)
 -- timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb.simpleioascii.filename <- ''
 -- Filename for spatial database not specified.
```

## Runtime Errors

```{code-block} pyrejournal
---
caption: Refer to [Step 1, Error 4](step01-error04.md).
---
 >> src/cig/pylith/libsrc/pylith/topology/MeshOps.cc:802:static void pylith::topology::MeshOps::checkMaterialLabels(const pylith::topology::Mesh&, pylith::int_array&)
 -- error (user-input)
 -- Material label_value '3' for cell '609' does not match the label_value of any materials or interfaces.
```

```{code-block} pyrejournal
---
caption: Refer to [Step 1, Error 5](step01-error05.md).
---
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Error converting spatial database values for gravitational_acceleration at (  -92285.3  -72772.7) in spatial database 'Gravity field'. Found near zero magnitude (0) for gravity field vector (  0  0).
```

```{code-block} pyrejournal
---
caption: Refer to [Step 6, Error 3](step06-error03.md).
---
 >> src/cig/pylith/libsrc/pylith/topology/FieldOps.cc:224:static void pylith::topology::FieldOps::checkSubfieldsExist(const pylith::string_vector &, const std::string &, const pylith::topology::Field &)
 -- error (user-input)
 -- Could not find 'lagrange_multiplier_fault' in domain solution field. Field contains: 'displacement'; the missing fields are required for interface 'FaultCohesiveKin'.
```

```{code-block} pyrejournal
---
caption: Refer to [Step 6, Error 4](step06-error04.md).
---
RuntimeError: Error occurred while reading spatial database file 'fault_slip.spatialdb'.
Read data for 3 out of 4 points.
```

```{code-block} pyrejournal
---
caption: Refer to [Step 6, Error 5](step06-error05.md).
---
RuntimeError: Could not find value 'final_slip_opening' in spatial database 'Fault rupture for main fault'. Available values are:
  final-slip-left-lateral
  final-slip-opening
  initiation-time
```

```{code-block} pyrejournal
---
caption: Refer to [Step 6, Error 6](step06-error06.md).
---
 >> src/cig/pylith/libsrc/pylith/topology/FieldQuery.cc:473:static void pylith::topology::_FieldQuery::findQueryIndices(FieldQuery::DBQueryContext *, const pylith::string_vector &)
 -- error (user-input)
 -- Could not find value 'final_slip_opening' in spatial database 'Fault rupture for main fault'. Available values are:
  final-slip-left-lateral
  final-slip-opening
  initiation-time
```

```{code-block} pyrejournal
---
caption: Refer to [Step 6, Error 7](step06-error07.md).
---
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Could not find values for initiation_time at (  -21860.3  -27621) in spatial database 'Fault rupture for main fault'.
```
