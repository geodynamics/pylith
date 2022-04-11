# PyLithApp

Full name: `pylith.apps.PyLithApp`

Python PyLithApp application.

## Pyre Facilities

* `dump_parameters`=\<component name\>: Dump parameters used and version information to file.
  - **current value**: 'dumpparamters', from {default}
  - **configurable as**: dumpparamters, dump_parameters
* `job`=\<component name\>: (no documentation available)
  - **current value**: 'job', from {default}
  - **configurable as**: job
* `launcher`=\<component name\>: (no documentation available)
  - **current value**: 'mpich', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/mpi/launchers/mpich.odb'} via {default}
  - **configurable as**: mpich, launcher
* `mesh_generator`=\<component name\>: Generates or imports the computational mesh.
  - **current value**: 'meshimporter', from {default}
  - **configurable as**: meshimporter, mesh_generator
* `metadata`=\<component name\>: Simulation metadata.
  - **current value**: 'metadata', from {default}
  - **configurable as**: metadata
* `petsc`=\<component name\>: Manager for PETSc options.
  - **current value**: 'petsc', from {default}
  - **configurable as**: petsc
* `problem`=\<component name\>: Computational problem to solve.
  - **current value**: 'timedependent', from {default}
  - **configurable as**: timedependent, problem
* `scheduler`=\<component name\>: (no documentation available)
  - **current value**: 'scheduler-none', from {file='/Users/baagaard/software/unix/py39-venv/pylith-debug/lib/python3.9/site-packages/pythia/pyre/schedulers/none.odb'} via {default}
  - **configurable as**: scheduler-none, scheduler
* `weaver`=\<component name\>: the pretty printer of my configuration as an XML document
  - **current value**: 'weaver', from {default}
  - **configurable as**: weaver

## Pyre Properties

* `include-citations`=\<bool\>: At end of simulation, display information on how to cite PyLith and components used.
  - **default value**: False
  - **current value**: False, from {default}
* `initialize_only`=\<bool\>: Stop simulation after initializing problem.
  - **default value**: False
  - **current value**: False, from {default}
* `nodes`=\<int\>: number of machine nodes
  - **default value**: 1
  - **current value**: 1, from {default}
* `start_python_debugger`=\<bool\>: Start python debugger at beginning of main().
  - **default value**: False
  - **current value**: False, from {default}
* `typos`=\<str\>: Specifies the handling of unknown properties and facilities
  - **default value**: 'pedantic'
  - **current value**: 'pedantic', from {default}
  - **validator**: (in ['relaxed', 'strict', 'pedantic'])
