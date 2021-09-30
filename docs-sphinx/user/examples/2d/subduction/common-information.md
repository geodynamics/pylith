# Common Information

As in previous examples, we place parameters common to the three steps in the `pylithapp.cfg` file so that we do not have to duplicate them for each step.
The settings contained in `pylithapp.cfg` for this problem consist of:

**pylithapp.journal.info** Settings that control the verbosity of the output written to stdout for the different components.

**pylithapp.mesh_generator** Settings that control mesh importing, such as the importer type, filename, and the spatial dimension of the mesh.

**pylithapp.timedependent** Settings that control the problem, such as the total time, time-step size, and spatial dimension.

**pylithapp.timedependent.materials** Settings that control the material type, specify which material IDs are to be associated with a particular material type, and give the name of the spatial database containing the physical properties for the material. The quadrature information is also given.

**pylithapp.problem.formulation.output** Settings related output of the solution over the domain and subdomain (ground surface).

**pylithapp.timedependent.materials.*MATERIAL*.output** Settings related to output of the state variables for material *MATERIAL*.

**pylithapp.petsc** PETSc settings to use for the problem such as the preconditioner type.

The physical properties for each material are specified in spatial database files.
For example, the elastic properties for the continental crust are in `mat_concrust.spatialdb`.
The provided spatial database files all use just a single point to specify uniform physical properties within each material.
A good exercise is to alter the spatial database files with the physical properties to match PREM.
