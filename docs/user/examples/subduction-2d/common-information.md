# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`, which contains numerous comments, so we only summarize the parameters here.

The settings contained in `pylithapp.cfg` for this problem consist of:

* `pylithapp.metadata` Metadata common to all of the simulations in the directory.
* `pylithapp.journal.info` Parameters that control the verbosity of the output written to stdout for the different components.
* `pylithapp.mesh_generator` Parameters for importing the finite-element mesh.
* `pylithapp.problem` Parameters that define the boundary value problem and it solution, such as the type of solver, solution fields.
* `pylithapp.problem.materials` Paramters that specify the governing equation and bulk rheologies.

The physical properties for each material are specified in spatial database files.
For example, the elastic properties for the continental crust are in `mat_concrust.spatialdb`.
The provided spatial database files all use just a single point to specify uniform physical properties within each material.
