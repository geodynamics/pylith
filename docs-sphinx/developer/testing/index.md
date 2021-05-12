# Testing

Debugging PyLith can be challenging due its many dependencies and complex interaction with PETSc.
In general, the most efficient strategy for debugging is to first try to expose the bug in a serial unit test, followed by an MMS Test, a serial full-scale test, and finally a parallel full-scale test.
This may require creating new tests if the bug is not exposed by current tests.
The PyLith developers make extensive use debuggers, such as `gdb` and `lldb`, and memory management analysis tools, such as `valgrind`, to detect and squash bugs. These are discussed in {ref}`developer-debugging-tools`.

```{toctree}
---
maxdepth: 2
---
libtests.md
mmstests.md
run-cxxtests.md
pytests.md
fullscale.md
debugging-tools.md
fields.md
ci-docker.md
```
