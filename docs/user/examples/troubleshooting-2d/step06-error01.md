# Step 6: Error 1

In troubleshooting Step 1 we resolved all of the errors associated with `pylithapp.cfg`.
In this example, we expect all errors to be associated with inputs files specific to Step 6, such as `step06_twofaults.cfg`.

## Error Message

```{code-block} pyrejournal
---
caption: Error message 1 when running Step 6.
linenos: True
emphasize-lines: 4-7
---
$ pylith step06_twofaults.cfg

# Output
 >> {default}::
 -- error (pyre.inventory)
 -- timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb.description <- ''
 -- Description for spatial database not specified.
 >> {default}::
 -- error (pyre.inventory)
 -- timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb.simpleioascii.filename <- ''
 -- Filename for spatial database not specified.
 >> step06_twofaults.cfg:81:
 -- error (pyre.inventory)
 -- timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb.filename <- 'fault_slip.spatialdb'
 -- unrecognized property 'timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb.filename'
usage: pylith [--<property>=<value>] [--<facility>.<property>=<value>] [FILE.cfg] ...
component 'pylithapp'
    properties: dump_parameters, help, help-components, help-persistence, help-properties, include-citations, initialize_only, job, launcher, metadata, nodes, petsc, problem, scheduler, show_application_flow, start_python_debugger, typos, weaver
    facilities: dump_parameters,job,launcher,metadata,petsc,problem,scheduler,weaver
For more information:
  --help-properties: prints details about user settable properties
  --help-components: prints details about user settable facilities and components
pylithapp: configuration error(s)
```

## Troubleshooting Strategy

We have several configuration errors, so we start with the first one at lines 3-6.
Spatial databases require a description, and PyLith cannot find one for `timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb`.
We examine the fault parameters for `step06_twofaults.cfg` and see that `db_auxiliary_field` for the earthquake rupture is missing the description.

## Resolution

```{code-block} cfg
---
caption: Correct error in `step06_twofaults.cfg`.
---
[pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
...
db_auxiliary_field.description = Fault rupture for main fault
```
