#!/bin/bash

# Run tests to generate coverage information. Upload test coverage data.
# Must run codecov script in top-level source directory.

make -j$(nproc) install

LCOV=`which lcov`
if test -f $LCOV; then
    make clean-coverage && \
	make -j$(nproc) check VERBOSE=1 -C tests/libtests; \
	if [ $? != 0 ]; then exit 1; fi; \
	make coverage-libtests && \
	pushd ../../src/pylith && \
	bash <(curl -s https://codecov.io/bash) -X gcov -f ../../build/pylith/coverage-libtests.info -F libtests -y ci-config/codecov.yml \
	    || echo "Codecov did not collect coverage reports" && popd
    
    #make coverage-pytests

    make clean-coverage && \
	make -j$(nproc) check VERBOSE=1 -C tests/mmstests; \
	if [ $? != 0 ]; then exit 1; fi; \
	make coverage-mmstests && \
	pushd ../../src/pylith && \
	bash <(curl -s https://codecov.io/bash) -X gcov -f ../../build/pylith/coverage-mmstests.info -F mmstests -y ci-config/codecov.yml \
	    || echo "Codecov did not collect coverage reports" && popd
    
    make clean-coverage && \
	make check VERBOSE=1 -C tests/fullscale; \
	if [ $? != 0 ]; then exit 1; fi; \
	make coverage-fullscale && \
	pushd ../../src/pylith && \
	bash <(curl -s https://codecov.io/bash) -X gcov -f ../../build/pylith/coverage-fullscale.info -F fullscale -y ci-config/codecov.yml \
	    || echo "Codecov did not collect coverage reports" && popd
    
else
    echo "lcov not found. Skipping test coverage."
fi
  
exit 0
