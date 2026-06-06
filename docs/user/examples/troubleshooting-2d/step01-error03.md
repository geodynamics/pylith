# Step 1: Error 3

## Error Message

```{code-block} pyrejournal
---
caption: Error message 3 when running Step 1.
linenos: True
emphasize-lines: 4-7
---
$ pylith step01a_gravity.cfg

# Output
 >> ./pylithapp.cfg:112:
 -- error (pyre.inventory)
 -- pylithapp.timedependent.materials.elasticity.auxiliary_subfields.bulk_modulus.basis_order <- '0'
 -- unknown component 'pylithapp.timedependent.materials.elasticity.auxiliary_subfields.bulk_modulus'
 >> ./pylithapp.cfg:113:
 -- error (pyre.inventory)
 -- pylithapp.timedependent.materials.elasticity.auxiliary_subfields.shear_modulus.basis_order <- '0'
 -- unknown component 'pylithapp.timedependent.materials.elasticity.auxiliary_subfields.shear_modulus'
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

Lines 3-6 indicate that we are trying to set the basis order for component `pylithapp.timedependent.materials.elasticity.auxiliary_subfields` that does not exist.
If we look up the documentation for components [`Elasticity`](../../components/materials/Elasticity.md) and [`IsotropicLinearElasticity`](../../components/materials/IsotropicLinearElasticity.md), then we see that the `bulk_modulus` auxiliary subfield is a subfield of the bulk rheology, not of the material.

## Resolution

```{code-block} cfg
---
caption: Correct error in `pylithapp.cfg`.
---
[pylithapp.problem.materials.wedge]
...
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0
```
