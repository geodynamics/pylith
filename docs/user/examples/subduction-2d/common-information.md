# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`, which contains numerous comments, so we only summarize the parameters here.

## Metadata, Mesh, and Output

The `pylithapp.metadata` section specifies metadata common to all simulations in the directory.
We control the verbosity of the output written to stdout using `journal.info`.
We set the parameters for importing the finite-element mesh in `pylithapp.mesh_generator`. 

## Physics

These quasi-static simulations solve the elasticity equation and include a fault, so we have a solution field with both displacement and Lagrange multiplier subfields.
%
\begin{gather}
\vec{s} = \left(\begin{array}{c} \vec{u} \quad \vec{\lambda} \end{array}\right)^T \\
\boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u}) = \vec{0}
\end{gather}
%
We use the default `TimeDependent` problem and solution field with a single displacement subfield of basis order 1.
In addition to output of the solution over the domain, we output the solution over the ground surface (+y boundary).
For the domain we skip one time step between writing the solution to the file.

```{code-block} cfg
---
caption: Solution and output parameters for all subduction-2d simulations.
---
[pylithapp.problem]
solution = pylith.problems.SolnDispLagrange

[pylithapp.problem]
solution_observers = [domain, groundsurf]
solution_observers.groundsurf = pylith.meshio.OutputSolnBoundary

[pylithapp.problem.solution_observers.domain]
# Skip 1 time step between output for the domain.
trigger.num_skip = 1

[pylithapp.problem.solution_observers.groundsurf]
# The `label` and `label_value` correspond to the name and tag of the
# physical group in the Gmsh Python script.
label = groundsurf
label_value = 10
```

The physical properties for each material are specified in spatial database files.
For example, the elastic properties for the continental crust are in `mat_concrust.spatialdb`.
The provided spatial database files all use just a single point to specify uniform physical properties within each material.

```{code-block} cfg
---
caption: Material parameters for the subduction-2d example suite. We only show the details for the continental crust material.
---
[pylithapp.problem]
materials = [continent_crust, ocean_crust, mantle]

[pylithapp.problem.materials]
[pylithapp.problem.materials.continent_crust]
description = Continental crust
label_value = 1

db_auxiliary_field.description = Continental crust properties
db_auxiliary_field.iohandler.filename = mat_concrust.spatialdb

observers.observer.trigger.num_skip = 1

auxiliary_subfields.density.basis_order = 0
bulk_rheology.auxiliary_subfields.bulk_modulus.basis_order = 0
bulk_rheology.auxiliary_subfields.shear_modulus.basis_order = 0
```
