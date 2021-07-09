// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "TestDriver.hh" // implementation of class methods

#include "pylith/topology/FieldOps.hh" // USES FieldOps::deallocate()
#include "pylith/utils/journals.hh" // USES journals.

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

namespace pylith {
    namespace testing {
        class _TestDriver {
public:

            /// Print help information.
            static
            void printHelp(void);

            /** Initialize PETSc.
             *
             * @param[in] argc Number of arguments passed.
             * @param[in] argv Array of input arguments.
             * @param[in] petscOptions Array of PETSc options to set.
             * @param[in] mallocDump Set malloc debug dump.
             */
            static
            int initializePetsc(int argc,
                                char* argv[],
                                const std::vector<std::string>& petscOptions,
                                const bool mallocDump);

            /** Add journal.
             *
             * @param[in] journals Array of journals.
             * @param[in] category Journal category.
             * @param[in] name Name of journal.
             */
            static
            void addJournal(TestDriver::journals_t& journals,
                            TestDriver::JournalEnum category,
                            const char* const name);

            /** Activate journals.
             *
             * @param[in] journals Journals to activate.
             */
            static
            void activateJournals(const TestDriver::journals_t& journals);

            /** List test hierarchy.
             *
             * @param[in] test Test to list.
             */
            static
            void printTests(const CppUnit::Test* const test);

            /** Find test matching name in test hierarchy.
             *
             * @param[in] test Test hierarchy.
             * @param[in] name Name of test to find.
             * @returns Test matching name or NULL if not found.
             */
            static
            const CppUnit::Test* findTest(const CppUnit::Test* const test,
                                          const std::string& name);

        };
    }
}

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
        _TestDriver::printHelp();
        return 0;
    } // if

    CppUnit::TestResultCollector result;
    try {
        // Initialize PETSc
        int err = _TestDriver::initializePetsc(argc, argv, _petscOptions, _mallocDump);CHKERRQ(err);

        // Initialize Python (to eliminate need to initialize when
        // parsing units in spatial databases).
        Py_Initialize();

        _TestDriver::activateJournals(_journals);

        CppUnit::Test* test = CppUnit::TestFactoryRegistry::getRegistry().makeTest();
        if (_listTests) {
            _TestDriver::printTests(test);
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
                const CppUnit::Test* testCase = _TestDriver::findTest(test, _tests[i]);
                if (testCase) {
                    runner.addTest(const_cast<CppUnit::Test*>(testCase));
                } else {
                    std::cerr << "ERROR: Could not find test '" << _tests[i] << "'." << std::endl;
                } // if
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
    static struct option options[9] = {
        {"help", no_argument, NULL, 'h'},
        {"list", no_argument, NULL, 'l'},
        {"quiet", no_argument, NULL, 'q'},
        {"tests", required_argument, NULL, 't'},
        {"petsc", required_argument, NULL, 'p'},
        {"journal.info", required_argument, NULL, 'i'},
        {"journal.debug", required_argument, NULL, 'd'},
        {"journal.warning", required_argument, NULL, 'w'},
        {0, 0, 0, 0}
    };

    _petscOptions.reserve(4);
    while (true) {
        // extern char* optarg;
        const char c = getopt_long(argc, argv, "hlqt:p:i:d:w:", options, NULL);
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
        case 'p': {
            if (_petscOptions.size()+1 > _petscOptions.capacity()) {
                _petscOptions.reserve(_petscOptions.capacity()+4);
            } // if
            _petscOptions.push_back(optarg);
            break;
        } // 'p'
        case 'i': {
            _TestDriver::addJournal(_journals, TestDriver::JOURNAL_INFO, optarg);
            break;
        } // 'p'
        case 'd': {
            _TestDriver::addJournal(_journals, TestDriver::JOURNAL_DEBUG, optarg);
            break;
        } // 'p'
        case 'w': {
            _TestDriver::addJournal(_journals, TestDriver::JOURNAL_WARNING, optarg);
            break;
        } // 'p'
        case '?':
            break;
        default:
            break;
        } // switch
    } // while
} // _parseArgs


// ---------------------------------------------------------------------------------------------------------------------
// Print help information.
void
pylith::testing::_TestDriver::printHelp(void) {
    std::cout << "Command line arguments:\n"
              << "[--help] [--list] [--quiet] [--tests=TEST_0,...,TEST_N\n\n"
              << "    --help            Print help information to stdout and exit.\n"
              << "    --list            Print names of tests.\n"
              << "    --quiet           Turn off dump of leaked memory.\n"
              << "    --tests           Comma separated list of tests to run (default is all tests).\n"
              << "    --petsc ARG=VALUE Arguments to pass to PETSc. May be repeated for multiple arguments.\n"
              << std::endl;
} // printHelp


// ---------------------------------------------------------------------------------------------------------------------
// Initialize PETSc.
int
pylith::testing::_TestDriver::initializePetsc(int argc,
                                              char* argv[],
                                              const std::vector<std::string>& petscOptions,
                                              const bool mallocDump) {
    int argcP = 1;
    char** argvP = new char*[1];
    argvP[0] = argv[0];
    PetscErrorCode err = PetscInitialize(&argcP, &argvP, NULL, NULL);CHKERRQ(err);
    delete[] argvP;argvP = NULL;

    if (mallocDump) {
        err = PetscOptionsSetValue(NULL, "-malloc_dump", "");CHKERRQ(err);
    } // if
    for (size_t i = 0; i < petscOptions.size(); ++i) {
        const size_t pos = petscOptions[i].find_first_of('=');
        if (pos < petscOptions[i].length()) {
            const std::string& arg = std::string("-") + petscOptions[i].substr(0, pos);
            const std::string& value = petscOptions[i].substr(pos+1);
            err = PetscOptionsSetValue(NULL, arg.c_str(), value.c_str());CHKERRQ(err);
        } else {
            const std::string& arg = std::string("-") + petscOptions[i];
            err = PetscOptionsSetValue(NULL, arg.c_str(), "");CHKERRQ(err);
        } // if/else
    } // for

    return 0;
} // initializePetsc


// ---------------------------------------------------------------------------------------------------------------------
// Add journal.
void
pylith::testing::_TestDriver::addJournal(TestDriver::journals_t& journals,
                                         TestDriver::JournalEnum category,
                                         const char* const name) {
    if (journals.size()+1 > journals.capacity()) {
        journals.reserve(journals.capacity()+4);
    } // if
    journals.push_back(std::pair<TestDriver::JournalEnum, std::string>(category, name));
} // addJournal


// ---------------------------------------------------------------------------------------------------------------------
// Activate journal.
void
pylith::testing::_TestDriver::activateJournals(const TestDriver::journals_t& journals) {
    for (size_t i = 0; i < journals.size(); ++i) {
        TestDriver::JournalEnum category = journals[i].first;
        const char* const name = journals[i].second.c_str();
        switch (category) {
        case TestDriver::JOURNAL_INFO: {
            pythia::journal::info_t journal(name);
            journal.activate();
            break;
        } // INFO
        case TestDriver::JOURNAL_DEBUG: {
            pythia::journal::debug_t journal(name);
            journal.activate();
            break;
        } // DEBUG
        case TestDriver::JOURNAL_WARNING: {
            pythia::journal::warning_t journal(name);
            journal.activate();
            break;
        } // INFO
        default:
            ;
            // PYLITH_JOURNAL_LOGICERROR("Unknown journal category '"<<category<<"'.");
        } // switch
    } // for
} // activateJournal


// ---------------------------------------------------------------------------------------------------------------------
/** List test hierarchy.
 *
 * @param[in] test Test to list.
 */
void
pylith::testing::_TestDriver::printTests(const CppUnit::Test* const test) {
    if (!test) { return; }
    std::cout << test->getName() << std::endl;
    if (!test->getChildTestCount()) { return; }
    for (int i = 0; i < test->getChildTestCount(); ++i) {
        _TestDriver::printTests(test->getChildTestAt(i));
    } // for
} // printTests


// ---------------------------------------------------------------------------------------------------------------------
// Find test matching name in test hierarchy.
const CppUnit::Test*
pylith::testing::_TestDriver::findTest(const CppUnit::Test* test,
                                       const std::string& name) {
    if (!test) { return NULL;}
    if (test->getName() == name) { return test; }
    if (!test->getChildTestCount()) { return NULL; }

    for (int i = 0; i < test->getChildTestCount(); ++i) {
        const CppUnit::Test* testCase = _TestDriver::findTest(test->getChildTestAt(i), name);
        if (testCase) { return testCase; }
    } // for

    return NULL;
} // findTest


// End of file
