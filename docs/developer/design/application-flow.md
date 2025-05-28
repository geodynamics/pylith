(sec-developer-design-application-flow)=
# PyLith Application Flow

The PyLith application driver performs two main functions.
First, it collects all user parameters from input files (e.g., `.cfg` files) and the command line, and then it performs from simple checks on the parameters.
Second, it launches the MPI job.

Once the MPI job launches, the application flow is:

1. Read the finite-element mesh; `pylith.meshio.MeshImporter`.
    1. Read the mesh (serial); `pylith::meshio::MeshIO`.
    2. Reorder the mesh, if desired; `pylith::topology::ReverseCuthillMcKee`.
    3. Insert cohesive cells as necessary (serial); `pylith::faults::FaultCohesive`.
    4. Distribute the mesh across processes (parallel); `pylith::topology::Distributor`.
    5. Refine the mesh, if desired (parallel); `pylith::topology::RefineUniform`.
2. Setup the problem.
    1. Preinitialize the problem by passing information from Python to C++ and doing minimal setup `pylith.Problem.preinitialize()`.
    2. Perform consistency checks and additional checks of user parameters; `pylith.Problem verifyConfiguration()`.
    3. Complete initialization of the problem; `pylith::problems::Problem::initialize()`.
3. Run the problem; `pylith.problems.Problem.run()`.
4. Cleanup; `pylith.problems.Problem.finalize()`.
    1. Close output files.
    2. Deallocate memory.
    3. Output PETSc log summary, if desired.

In the first step, we list the object performing the work, whereas in subsequent steps we list the top-level object method responsible for the work.
Python objects are listed using the `path.class` syntax while C++ objects are listed using `namespace::class` syntax.
Note that a child class may redefine or perform additional work compared to what is listed in the parent class method.

Reading the mesh and the first two steps of the problem setup are controlled from Python.
That is, at each step Python calls the corresponding C++ methods using SWIG.
Starting with the complete initialization of the problem, the flow is controlled at the C++ level.

## Time-Dependent Problem

In a time-dependent problem the PETSc `TS` object (relabeled `PetscTS` within PyLith) controls the time stepping.
Within each time step, the `PetscTS` object calls the PETSc linear and nonlinear solvers as needed, which call the following methods of the C++ `pylith::problems::TimeDependent` object as needed: `computeRHSResidual()`, `computeLHSResidual()`, and `computeLHSJacobian()`.
The `pylith::problems::TimeDependent` object calls the corresponding methods in the boundary conditions, constraints, and materials objects.
At the end of each time step, it calls `problems::TimeDependent::poststep()`.  
