# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`, which contains numerous comments, so we only summarize the parameters here.

The settings contained in `pylithapp.cfg` for this problem consist of:

:`pylithapp.metadata`: Metadata common to all of the simulations in the directory.
:`pylithapp.mesh_generator`: Parameters for importing the finite-element mesh.
:`pylithapp.problem`: Parameters that define the boundary value problem and its solution, such as the type of solver, solution fields, and output over the domain.
:`pylithapp.problem.materials`: Parameters that specify the governing equation and bulk rheologies.
:`pylithapp.problem.bc`: Parameters that specify the boundary conditions.

To make it easier to switch between different bulk rheologies for different simulations, we separate the material parameters relevant to the bulk rheologies into three parameter files: `mat_elastic.cfg`, `mat_viscoelastic.cfg`, and `mat_elastic_incompressible.cfg`.
We use `mat_elastic.cfg` when we want all materials to use elastic bulk rheologies.
We use `mat_viscoelastic.cfg` when we want the mantle and bottom of the slab to use viscoelastic bulk rheologies. We use `mat_elastic_incompressible.cfg` only for Step 8b in which we use incompressible elasticity.
The physical properties for each material are specified in spatial database files.
For example, the elastic properties for the crust are in `mat_crust_elastic.spatialdb`.
The spatial database files for elastic properties all use just a single point to specify uniform physical properties within each material.

The distorted cells at the trench result in a small degradation in the preconditioning, so the default KSP absolute tolerance is not quite small enough for the nonlinear solver in a single iteration.
Consequently, we reduce the KSP absolute tolerance to 2.0e-8, which allow the nonlinear solver to converge in a single iteration.
