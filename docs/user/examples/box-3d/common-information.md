# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`.
This limits duplicate information in the parameter files for each simulation.
The `pylithapp.cfg` file contains numerous comments, so we only summarize the parameters here.

The settings contained in `pylithapp.cfg` for this problem consist of:

* `pylithapp.metadata` Metadata common to all of the simulations in the directory.
* `pylithapp.journal.info` Parameters that control the verbosity of the output written to stdout for the different components.
* `pylithapp.mesh_generator` Parameters for importing the finite-element mesh.
* `pylithapp.problem` Parameters that define the boundary value problem and it solution, such as the type of solver, solution fields.
* `pylithapp.problem.materials` Paramters that specify the governing equation and bulk rheologies.

The physical properties for each material are specified in a spatial database.
In this case we have only one material and the material properties are uniform, so we use a single `UniformDB` to specify the material properties in `pylithapp.cfg`.
When we have a spatial variation in the material properties, we specify them in a spatial database file.
