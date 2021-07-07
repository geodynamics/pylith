# Testing

Testing and debugging PyLith can be challenging due its many dependencies and complex interaction with PETSc.
Our strategy is to test at a variety of levels to isolate bugs close to their origin while also building a comprehensive suite of tests.
We use C++ unit tests for verifying small pieces of code, for example individual C++ class methods,  work as intended.
We use the Method of Manufactured Solutions for verifying implementation of governing equations and their solution using PETSc.
We use full-scale tests for verifying integration for complete simulations in both serial and parallel.
When these tests, an example, or user simulation suggests a bug exists, we leverage the test suite.
In general, the most efficient strategy for debugging is to first try to expose the bug in a serial unit test, followed by an MMS Test, a serial full-scale test, and finally a parallel full-scale test.
This may require creating new tests if the bug is not exposed by current tests.
The PyLith developers make extensive use of debuggers, such as `gdb` and `lldb`, and memory management analysis tools, such as `valgrind`, to detect and squash bugs. These are discussed in {ref}`developer-debugging-tools`.

```{toctree}
libtests.md
mmstests.md
run-cxxtests.md
pytests.md
fullscale.md
debugging-tools.md
fields.md
ci-docker.md
```
