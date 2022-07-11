# Step 1: Error 2

## Error Message

```{code-block} console
---
caption: Error message 2 when running Step 1.
linenos: True
emphasize-lines: 3-6
---
$ pylith step01_gravity.cfg

 >> {default}::
 -- pyre.inventory(error)
 -- timedependent.problem_defaults.name <- ''
 -- Missing required property 'name' in default options for problem.
 >> ./pylithapp.cfg:133:
 -- pyre.inventory(error)
 -- pylithapp.timedependent.materials.elasticity.auxiliary_subfields.bulk_modulus.basis_order <- '0'
 -- unknown component 'pylithapp.timedependent.materials.elasticity.auxiliary_subfields.bulk_modulus'
 >> ./pylithapp.cfg:134:
 -- pyre.inventory(error)
 -- pylithapp.timedependent.materials.elasticity.auxiliary_subfields.shear_modulus.basis_order <- '0'
 -- unknown component 'pylithapp.timedependent.materials.elasticity.auxiliary_subfields.shear_modulus'
usage: pylith [--<property>=<value>] [--<facility>.<property>=<value>] [FILE.cfg] ...
component 'pylithapp'
    properties: dump_parameters, help, help-components, help-persistence, help-properties, include-citations, initialize_only, job, launcher, mesh_generator, metadata, nodes, petsc, problem, scheduler, start_python_debugger, typos, weaver
    facilities: dump_parameters,job,launcher,mesh_generator,metadata,petsc,problem,scheduler,weaver
For more information:
  --help-properties: prints details about user settable properties
  --help-components: prints details about user settable facilities and components
pylithapp: configuration error(s)
```

## Troubleshooting Strategy

The full error message shows several errors associated with our parameter files.
Because errors often cascade, the best place to start is to address the first error (lines 3-6).
We see that the default name for `problem_defaults` is empty and PyLith requires one.
This name will be used as the first part of the filename for output.
We resolve this error by examining the parameters for `pylithapp.timedependent.problem_defaults` in `step01_gravity.cfg`.
We see that we do not set any of the problem defaults.

## Resolution

We must set the simulation name in the problem defaults.

```{code-block} cfg
---
caption: Correct error in `pylithapp.cfg`.
---
[pylithapp]
...
# Set the name of the problem that will be used to construct the
# output filenames. The default directory for output is 'output'.
problem.defaults.name = step01_gravity
```
