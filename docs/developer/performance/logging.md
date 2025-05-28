# Logging

PyLith includes `pylith::utils::EventLogger` as a high-level interface for event logging with PETSc.
We create local `Events` classes within the C++ implementation files (`*.cc`) to create variables that hold the event ids.
This streamlines logging events within the C++ code.
We identify events using variables rather than strings, which means that typos in names are detected at compile time rather than at runtime.
Refer to the `_Problem::Events` class near the top of `libsrc/pylith/problems/Problem.cc` for an example of how we implement event logging.

We name events in PyLith C++ code following the template `PL:CLASS:METHOD` (C++ namespace style but reducing `::` to `:` to reduce the length of the string) in which we are logging method `METHOD` in C++ class `CLASS`.
In the PyLith Python code, we follow the template `PL.CLASS.METHOD`.
