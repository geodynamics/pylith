# Running tests in CI Docker containers

If one or more of the CI test runners report an error or someone reports an error using a specific Linux distribution, it is often convenient to run a CI test runner interactively on your local machine using Docker. Currently, we maintain Docker images containing all of the PyLith dependencies for the following Linux distributions:

* debian-stable
* debian-testing
* ubuntu-18.04
* ubuntu-20.04
* ubuntu-20.10
* fedora-32
* fedora-33
* centos-7
* centos-8
  
You need to checkout the PyLith branch that you want to test and be in the top-level PyLith source directory.

```{code-block} console
# Change to the top-level PyLith source directory.
cd PYLITH_TOP_SRC_DIR

# Set the name of the Linux distribution to use to the BASE_IMAGE environment
# variable.
$ BASE_IMAGE=debian-stable

# Build a new Docker image with the PyLith source code. The base image will
# be downloaded as necessary. The "--target build" command line argument
# will run "make install" but not "make check". To stop at the configure step,
# use "--target src".
$ docker build \
  -t pylith-debug \
  --build-arg BASE_IMAGE=registry.gitlab.com/cig-pylith/pylith_installer/test-env/$BASE_IMAGE \
  --target build \
  -f docker/pylith-testenv .

# Run the Docker image interactively. This allows you to run "make check" or
# run tests manually, potentially making use of the gdb debugger and valgrind.
$ docker run -ti --rm pylith-debug /bin/bash
