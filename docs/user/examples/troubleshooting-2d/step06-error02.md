# Step 6: Error 2

## Error Message

```{code-block} pyrejournal
---
caption: Error message 2 when running Step 6.
linenos: True
emphasize-lines: 4-7
---
$ pylith step06_twofaults.cfg

# Output
 >> {default}::
 -- error (pyre.inventory)
 -- timedependent.interfaces.faultcohesivekin.singlerupture.kinsrcstep.simpledb.simpleioascii.filename <- ''
 -- Filename for spatial database not specified.
 >> step06_twofaults.cfg:84:
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

We still have problem configuration errors.
The filename for the `SimpleIOAscii` reader for a `SimpleDB` is missing.
We examine the earthquake rupture parameters in `step06_twofaults.cfg` again and see that the filename is associated with the `SimpleDB`; however, if we look at the [`SimpleDB` documentation](https://spatialdata.readthedocs.io/en/latest/user/components/spatialdb/SimpleDB.html) we notice that it has a reader `SimpleIOAscii`, which has a filename.

## Resolution

```{code-block} cfg
---
caption: Correct error in `step06_twofaults.cfg`.
---
# Error
db_auxiliary_field.filename = fault_slip.spatialdb

# Correct
db_auxiliary_field.iohandler.filename = fault_slip.spatialdb
```
