(sec-developer-design-code-layout)=
# Code Layout

## Directory Structure

The C++, Python, and SWIG Python/C++ interface files all sit in different directories.
Similarly, the unit tests, MMS tests, full-scale tests, examples and documentation are also in their own directories.

```bash
pylith/
├── ci-config # continuous integration testing configuration
├── docker # Dockerfiles for creating Docker images
├── developer # Utilities for developers
├── docs # Documentation source code
├── libsrc # C++ source for library
├── modulesrc # SWIG interfaces files for C++ Python bindings.
├── pylith # Python source code for PyLith modules.
├── applications # Source code for command line programs
├── m4 # Autoconf macros
├── share # Common parameter settings
├── examples # Example suite
└── tests
    ├── libtests # C++ unit tests
    ├── mmstests # C++ Method of Manufactured solution tests
    ├── pytests # Python unit tests
    ├── fullscale # Automated full-scale tests
    └── manual # Manual full-scale tests
```

The C++ library, SWIG interface files, and Python modules are organized into several subpackages with common names.

* **bc**  Boundary conditions.
* **faults** Faults.
* **feassemble** General finite-element formulation.
* **fekernels** Finite-element pointwise functions (kernels).
* **friction** Fault constitutive models.
* **materials** Material behavior, including bulk constitutive models.
* **meshio** Input and output.
* **problems** General problem formulation.
* **testing** Common testing infrastructure.
* **topology** Finite-element mesh topology.
* **utils** General utilities.
