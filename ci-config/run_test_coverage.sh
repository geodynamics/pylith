#!/bin/bash

# Run tests to generate coverage information. Upload test coverage data.
# Must run codecov script in top-level source directory.

# ------------------------------------------------------------------------------
# Build
# ------------------------------------------------------------------------------
make -j$(nproc) install
if [ $? != 0 ]; then exit 1; fi


# ------------------------------------------------------------------------------
# libtests
# ------------------------------------------------------------------------------
make clean-coverage && make -j$(nproc) check -C tests/libtests VERBOSE=1
if [ $? != 0 ]; then exit 1; fi

make coverage-libtests
if [ $? != 0 ]; then exit 1; fi

if [ -r coverage-libtests.info ]; then
  pushd ../../src/spatialdata && \
      bash <(curl -s https://codecov.io/bash) -X gcov -f ../../build/spatialdata/coverage-libtests.info -F libtests -y ci-config/codecov.yml \
	  || echo "Codecov did not collect coverage reports." && \
      popd
fi


# ------------------------------------------------------------------------------
# mmstests
# ------------------------------------------------------------------------------
make clean-coverage && make -j$(nproc) check -C tests/mmstests VERBOSE=1
if [ $? != 0 ]; then exit 1; fi

make coverage-mmstests
if [ $? != 0 ]; then exit 1; fi

if [ -r coverage-mmstests.info ]; then
  pushd ../../src/spatialdata && \
      bash <(curl -s https://codecov.io/bash) -X gcov -f ../../build/spatialdata/coverage-mmstests.info -F mmstests -y ci-config/codecov.yml \
	  || echo "Codecov did not collect coverage reports." && \
      popd
fi


# ------------------------------------------------------------------------------
# pytests
# ------------------------------------------------------------------------------
#make -j$(nproc) check -C tests/pytests VERBOSE=1
#if [ $? != 0 ]; then exit 1; fi

#make coverage-pytests
#if [ $? != 0 ]; then exit 1; fi

#if [ -r coverage-pytests.xml ]; then
#  pushd ../../src/spatialdata && \
#      bash <(curl -s https://codecov.io/bash) -X gcov -f ../../build/spatialdata/coverage-pytests.xml -F pytests -y ci-config/codecov.yml \
#	  || echo "Codecov did not collect coverage reports." && \
#      popd
#fi


# ------------------------------------------------------------------------------
# fullscale
# ------------------------------------------------------------------------------
make clean-coverage && make -j$(nproc) check -C tests/fullscale VERBOSE=1
if [ $? != 0 ]; then exit 1; fi

make coverage-fullscale
if [ $? != 0 ]; then exit 1; fi

if [ -r coverage-fullscale.info ]; then
  pushd ../../src/spatialdata && \
      bash <(curl -s https://codecov.io/bash) -X gcov -f ../../build/spatialdata/coverage-fullscale.info -F fullscale -y ci-config/codecov.yml \
	  || echo "Codecov did not collect coverage reports." && \
      popd
fi


exit 0
