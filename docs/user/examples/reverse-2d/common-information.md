# Common Information

In addition to the finite-element mesh, PyLith requires files to specify the simulation parameters.
We specify parameters common to all simulations in a directory in `pylithapp.cfg`, which contains numerous comments, so we only summarize the parameters here.

The settings contained in `pylithapp.cfg` for this problem consist of:

* `pylithapp.metadata` Metadata common to all of the simulations in the directory.
* `pylithapp.journal.info` Parameters that control the verbosity of the output written to stdout for the different components.
* `pylithapp.mesh_generator` Parameters for importing the finite-element mesh.
* `pylithapp.problem` Parameters that define the boundary value problem and it solution, such as the type of solver, solution fields.
* `pylithapp.problem.materials` Paramters that specify the governing equation and bulk rheologies.
* `pylithapp.problem.bc` Boundary condition parameters common aomong the simulations in this directory.

These static and quasistatic simulations solve the elasticity equation.
We use the default solution field (displacement); we will override this for the simulations with a fault.

We use the same material properties for several simulations in this directory, so we specify them in `pylithapp.cfg` to avoid repeating the information in the file with parameters for each simulation.
We use a `SimpleDB` spatial database so that we can simply use a different spatial database file when we change the bulk rheology.
Similarly, for all of the simulations in this directory we use Dirichlet (displacement) boundary conditions on the +x, -x, and -y boundaries that constrain the displacement component perpendicular to the fault.
