# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`, which contains numerous comments, so we only summarize the parameters here.

The settings contained in `pylithapp.cfg` for this problem consist of:

* `pylithapp.metadata` Metadata common to all of the simulations in the directory.
* `pylithapp.journal.info` Parameters that control the verbosity of the output written to stdout for the different components.
* `pylithapp.mesh_generator` Parameters for importing the finite-element mesh.
* `pylithapp.problem` Parameters that define the boundary value problem and it solution, such as the type of solver, solution fields, and output over the domain.
* `pylithapp.problem.materials` Parameters that specify the governing equation and bulk rheologies.
* `pylithapp.problem.bc` Parameters that specify the boundary conditions.

To make it easier to switch between elastic and viscoelastic bulk rheologies for different simulations, we separate the material parameters relevant to the bulk rheologies into two parameter files: `mat_elastic.cfg` and `mat_viscoelastic.cfg`.
We use `mat_elastic.cfg` when we want all materials to use elastic bulk rheologies.
We use `mat_viscoelastic.cfg` when we want the mantle and bottom of the slab to use viscoelastic bulk rheologies.
The physical properties for each material are specified in spatial database files.
For example, the elastic properties for the crust are in `mat_crust_elastic.spatialdb`.
The spatial database files for elastic properties all use just a single point to specify uniform physical properties within each material.
