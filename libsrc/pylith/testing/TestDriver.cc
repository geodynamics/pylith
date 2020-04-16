// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "TestDriver.hh" // implementation of class methods

#include "pylith/topology/FieldOps.hh" // USES FieldOps::deallocate()

#include "petsc.h"
#include <Python.h>

#include <cppunit/extensions/TestFactoryRegistry.h>

#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/TextOutputter.h>

#include <stdlib.h> // USES abort()
#include <getopt.h> // USES getopt_long()
#include <sstream> // USES std::ostringstream, std::istringstream

// ----------------------------------------------------------------------
// Constructor
pylith::testing::TestDriver::TestDriver() :
    _showHelp(false),
    _listTests(false),
    _mallocDump(true) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::testing::TestDriver::~TestDriver(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Run info application.
int
pylith::testing::TestDriver::run(int argc,
                                 char* argv[]) {
    _parseArgs(argc, argv);

    if (_showHelp) {
        _printHelp();
        return 0;
    } // if

    CppUnit::TestResultCollector result;
    try {
        // Initialize PETSc
        PetscErrorCode err = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(err);
        if (_mallocDump) {
            err = PetscOptionsSetValue(NULL, "-malloc_dump", "");CHKERRQ(err);
        } // if

        // Initialize Python (to eliminate need to initialize when
        // parsing units in spatial databases).
        Py_Initialize();

        CppUnit::Test* test = CppUnit::TestFactoryRegistry::getRegistry().makeTest();
        if (_listTests) {
            _printTests(test);
            return 0;
        } // if

        CppUnit::TestResult controller;
        CppUnit::BriefTestProgressListener progress;
        controller.addListener(&result);
        controller.addListener(&progress);

        CppUnit::TestRunner runner;
        if (_tests.empty()) {
            runner.addTest(test);
        } else {
            for (size_t i = 0; i < _tests.size(); ++i) {
                const CppUnit::Test* testCase = _findTest(test, _tests[i]);
                runner.addTest(const_cast<CppUnit::Test*>(testCase));
            } // for
        } // if/else
        runner.run(controller);

        // Print test results
        CppUnit::TextOutputter outputter(&result, std::cerr);
        outputter.write();

        // Finalize Python
        Py_Finalize();

        // Finalize PETSc
	pylith::topology::FieldOps::deallocate();
        err = PetscFinalize();CHKERRQ(err);
    } catch (...) {
        abort();
    } // catch

    if (!_mallocDump) {
        std::cout << "WARNING -malloc dump is OFF\n" << std::endl;
    } // if

    return (result.wasSuccessful() ? 0 : 1);
} // run


// ---------------------------------------------------------------------------------------------------------------------
// Parse command line arguments.
void
pylith::testing::TestDriver::_parseArgs(int argc,
                                        char* argv[]) {
    static struct option options[5] = {
        {"help", no_argument, NULL, 'h'},
        {"list", no_argument, NULL, 'l'},
        {"quiet", no_argument, NULL, 'q'},
        {"tests", required_argument, NULL, 't'},
        {0, 0, 0, 0}
    };

    while (true) {
        // extern char* optarg;
        const char c = getopt_long(argc, argv, "hlm:", options, NULL);
        if (-1 == c) { break; }
        switch (c) {
        case 'h':
            _showHelp = true;
            break;
        case 'l':
            _listTests = true;
            break;
        case 'q':
            _mallocDump = false;
            break;
        case 't': {
            _tests.clear();
            std::istringstream tokenStream(optarg);
            std::string token;
            while (std::getline(tokenStream, token, ',')) {
                _tests.push_back(token);
            } // while
            break;
        } // 't'
        case '?':
            break;
        default:
            // May have PETSc command line arguments.
            break;
        } // switch
    } // while
} // _parseArgs


// ---------------------------------------------------------------------------------------------------------------------
// Print help information.
void
pylith::testing::TestDriver::_printHelp(void) {
    std::cout << "Command line arguments:\n"
              << "[--help] [--list] [--quiet] [--tests=TEST_0,...,TEST_N\n\n"
              << "    --help           Print help information to stdout and exit.\n"
              << "    --list           Print names of tests.\n"
              << "    --quiet          Turn off dump of leaked memory.\n"
              << "    --tests          Comma separated list of tests to run (default is all tests).\n"
              << std::endl;
} // _printHelp


// ---------------------------------------------------------------------------------------------------------------------
/** List test hierarchy.
 *
 * @param[in] test Test to list.
 */
void
pylith::testing::TestDriver::_printTests(const CppUnit::Test* const test) {
    if (!test) { return; }
    std::cout << test->getName() << std::endl;
    if (!test->getChildTestCount()) { return; }
    for (int i = 0; i < test->getChildTestCount(); ++i) {
        _printTests(test->getChildTestAt(i));
    } // for
} // _printTests


// ---------------------------------------------------------------------------------------------------------------------
/** Find test matching name in test hierarchy.
 *
 * @param[in] test Test hierarchy.
 * @param[in] name Name of test to find.
 * @returns Test matching name or NULL if not found.
 */
const CppUnit::Test*
pylith::testing::TestDriver::_findTest(const CppUnit::Test* test,
                                       const std::string& name) {
    if (!test) { return NULL;}
    if (test->getName() == name) { return test; }
    if (!test->getChildTestCount()) { return NULL; }

    for (int i = 0; i < test->getChildTestCount(); ++i) {
        const CppUnit::Test* found = _findTest(test->getChildTestAt(i), name);
        if (found) { return found; }
    } // for

    return NULL;
} // _findTest


// End of file
