# Rebuilding PETSc and PyLith

For instructions on how to build PyLith from source, please see the [PyLith Installer Documentation](https://pylith-installer.readthedocs.io/en/latest/devenv/index.html).

## Rebuilding PETSc

Updating and rebuilding PETSc is quite simple once it has been configured and built once before.

```{code-block} bash
# Change to PETSc source Directory
cd $PETSC_DIR

# Fetch the updates, pruning deleted branches.
git fetch -p

# Ensure you are on the `knepley/pylith` branch and then pull updates.
git checkout knepley/pylith
git pull

# Reconfigure
arch-pylith-debug/lib/petsc/conf/reconfigure-arch-pylith-debug.py

# Rebuild
make
```

:::{important}
After rebuilding PETSc, you should rebuild PyLith.
If there are incompatibilities between the two, then you will normally get compilation errors building PyLith.
If there are incompatibilities and you do not rebuild PyLith, then you will usually get a segmentation fault when running PyLith.
:::

## Rebuilding PyLith

### Overview

Pylith uses the GNU Build System (often called autotools, which consists of autoconf, automake, and libtool) to configure and build.
The configure options and checks are defined in `configure.ac` with additional macros in the `m4` directory.
Note that the `m4` directory is a Git submodule corresponding to the `geodynamics/autoconf_cig` Git repository.

### Updating your fork

As long as you setup your `main` branch to correspond to `upstream/main`, you should never have to use the `main` branch in your fork.
You only use your fork for your feature branches.

### Makefiles

The PyLith `configure` script uses automake to convert each `Makefile.am` file into a `Makefile`.
The organization and content of the `Makefile.am` file depends on whether it is related to the C++ library, SWIG interface files, Python modules, C++ unit tests, Python unit tests, MMS tests, full-scale tests, or examples.

For the C++ library files within the `libsrc` directory, the `libsrc/Makefile.am` contains the implementation files while the header files are listed in the `Makefile.am` file within the underlying directories.
For the SWIG interface files within the `modulesrc` directory, `Makefile.am` file contains information on how to build the SWIG Python module.
The Python modules use a single `Makefile.am` file in the `pylith` directory.
Each test suite (C++ unit tests for each subpackage, Python unit tests, and each suite of full-scale tests) use a single `Makefile.am`.
These files define how the tests are built, additional input files that should be included in the source distribution, and temporary files that should be deleted.

Each suite of examples contains a `Makefile.am` that defines the files that are to be included in the source distribution.
It also defines which files are created and should be deleted upon `make clean`.

### Build Targets

Several build targets are defined to make it easy to rebuild/reinstall PyLith and rerun tests whenever the source code changes.
Each target can be run using `make TARGET` where `TARGET` is one of the following:

`all`
: Build all source code.

`install`
: Build and install.

`check`
: Run the entire test suite.

:::{tip}
On a machine with multiple cores, faster builds of the C++ code (library and C++ unit tests) are available using the `-jNTHREADS` command line argument, where `NTHREADS` is the number of threads to use.
We usually set the number of threads equal to twice the number of physical cores.
:::

After modifying code, the C++ library, SWIG modules, and Python code need to be rebuilt and reinstalled before running a PyLith simulation.

```{code-block} bash
---
caption: Rebuilding Pylith C++ library, SWIG modules, and Python modules
---
# Change to top-level PyLith build directory.
cd $PYLITH_BUILDDIR

# Reinstall everything using 16 threads to build library.
make install -j16

# Rebuild and reinstall only the library using 16 threads.
make install -j16 -C libsrc

# Rebuild and reinstall only the SWIG modules
make install -C modulesrc

# Reinstall only the Python modules
make install -C pylith
```

After modifying the C++ code, only the C++ library needs to be rebuilt before running C++ unit tests or MMS tests.

```{code-block} bash
---
caption: Rebuilding Pylith C++ library and rerunning the C++ unit tests
---
# Change to top-level PyLith build directory.
cd $PYLITH_BUILDDIR/build/pylith-debug

# Rebuild the library using 16 threads.
make -j16 -C libsrc

# Rerun the C++ unit tests.
make check -C tests/libtests
```

Similarly, after changing the Python code, only the Python modules need to be reinstalled before running a Python unit test.
However, if changes are also made to the C++ code, then the C++ library and modules must be rebuilt and reinstalled before running a Python unit test.

```{code-block} bash
---
caption: Rebuilding Pylith Python modules and rerunning Python unit tests
---
# Change to top-level PyLith build directory.
cd $PYLITH_BUILDDIR/build/pylith-debug

# Reinstall PyLith Python modules.
make install -C pylith

# Rerun Python unit tests.
make check -C tests/pytests
```
