# Debugging

Debugging PyLith can be challenging due its many dependencies and complex interaction with PETSc.
In general, the most efficient strategy for debugging is to first try to expose the bug in a serial unit test, followed by an MMS Test, a serial full-scale test, and finally a parallel full-scale test.
This may require creating new tests if the bug is not exposed by current tests.
The PyLith developers make extensive use debuggers, such as `gdb` and `lldb`, and memory management analysis tools, such as `valgrind`, to detect and squash bugs.

:::{toctree}
---
maxdepth: 2
---
debugger.md
valgrind.md
fields.md
:::
