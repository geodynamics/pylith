// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/testing/testingfwd.hh"
#include "pylith/topology/FieldOps.hh" // USES FieldOps::deallocate()
#include "pylith/utils/journals.hh" // USES journals

#include "petsc.h"
#include <Python.h>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/Test.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/TextOutputter.h>

#include <getopt.h> // USES getopt_long()
#include <cstdlib> // USES abort()
#include <sstream> // USES std::ostringstream, std::istringstream
#include <vector> // USES std::vector
#include <utility> // USES std::pair
#include <string> // USES std::string

// ------------------------------------------------------------------------------------------------
class pylith::testing::TestDriver {
public:

    enum JournalEnum {
        JOURNAL_INFO=0,
        JOURNAL_DEBUG=1,
        JOURNAL_WARNING=2,
    };

    typedef std::vector< std::pair<JournalEnum,std::string> > journals_t;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDriver(void);

    /** Run test application.
     * @param argc[in] Number of arguments passed.
     * @param argv[in] Array of input arguments.
     *
     * @returns 1 if errors were detected, 0 otherwise.
     */
    int run(int argc,
            char* argv[]);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /** Parse command line arguments.
     *
     * @param[in] argc Number of arguments passed.
     * @param[in] argv Array of input arguments.
     */
    void _parseArgs(int argc,
                    char* argv[]);

    /// Print help information.
    void _printHelp(void);

    /** List test hierarchy.
     *
     * @param[in] test Test to list.
     */
    void _printTests(const CppUnit::Test* const test);

    /** Initialize PETSc.
     *
     * @param[in] argc Number of arguments passed.
     * @param[in] argv Array of input arguments.
     * @param[in] petscOptions Array of PETSc options to set.
     * @param[in] mallocDump Set malloc debug dump.
     */
    int _initializePetsc(int argc,
                         char* argv[],
                         const std::vector<std::string>& petscOptions,
                         const bool mallocDump);

    /** Add journal.
     *
     * @param[in] category Journal category.
     * @param[in] name Name of journal.
     */
    void _addJournal(JournalEnum category,
                     const char* const name);

    /// Activate journals.
    void _activateJournals(void);

    /** Find test matching name in test hierarchy.
     *
     * @param[in] test Test hierarchy.
     * @param[in] name Name of test to find.
     * @returns Test matching name or NULL if not found.
     */
    const CppUnit::Test* _findTest(const CppUnit::Test* const test,
                                   const std::string& name);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::vector<std::string> _tests;
    std::vector<std::string> _petscOptions;
    journals_t _journals;
    bool _showHelp;
    bool _listTests;
    bool _mallocDump;

    // NOT IMPLEMENTED
    // ////////////////////////////////////////////////////////////////////////////////////////////
private:

    TestDriver(const TestDriver&); ///< Not implemented
    const TestDriver& operator=(const TestDriver&); ///< Not implemented

};

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::testing::TestDriver::TestDriver(void) :
    _showHelp(false),
    _listTests(false),
    _mallocDump(true) {}


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
        int err = _initializePetsc(argc, argv, _petscOptions, _mallocDump);CHKERRQ(err);

        // Initialize Python (to eliminate need to initialize when
        // parsing units in spatial databases).
        Py_Initialize();

        CppUnit::Test* test = CppUnit::TestFactoryRegistry::getRegistry().makeTest();
        if (_listTests) {
            _printTests(test);
            return 0;
        } // if

        _activateJournals();

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
        std::cout << "WARNING PETSc option -malloc dump is OFF\n" << std::endl;
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
        const int c = getopt_long(argc, argv, "hlqt:p:i:d:w:", options, NULL);
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
            _addJournal(JOURNAL_INFO, optarg);
            break;
        } // 'p'
        case 'd': {
            _addJournal(JOURNAL_DEBUG, optarg);
            break;
        } // 'p'
        case 'w': {
            _addJournal(JOURNAL_WARNING, optarg);
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
pylith::testing::TestDriver::_printHelp(void) {
    std::cout << "Command line arguments:\n"
              << "[--help] [--list] [--quiet] [--tests=TEST_0,...,TEST_N\n\n"
              << "    --help            Print help information to stdout and exit.\n"
              << "    --list            Print names of tests.\n"
              << "    --quiet           Turn off dump of leaked memory.\n"
              << "    --tests           Comma separated list of tests to run (default is all tests).\n"
              << "    --petsc ARG=VALUE Arguments to pass to PETSc. May be repeated for multiple arguments.\n"
              << "    --journal.TYPE=COMPONENT Activate journal (TYPE=[info,debug,warning,error]) for COMPONENT. May be repeated with different components.\n"
              << std::endl;
} // _printHelp


// ---------------------------------------------------------------------------------------------------------------------
// Initialize PETSc.
int
pylith::testing::TestDriver::_initializePetsc(int argc,
                                              char* argv[],
                                              const std::vector<std::string>& petscOptions,
                                              const bool mallocDump) {
    int argcP = 1;
    char** argvP = new char*[argcP+1];
    argvP[0] = argv[0];
    argvP[argcP] = NULL; // C standard is argv[argc] == NULL.
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
} // _initializePetsc


// ---------------------------------------------------------------------------------------------------------------------
// Add journal.
void
pylith::testing::TestDriver::_addJournal(JournalEnum category,
                                         const char* const name) {
    if (_journals.size()+1 > _journals.capacity()) {
        _journals.reserve(_journals.capacity()+4);
    } // if
    _journals.push_back(std::pair<JournalEnum, std::string>(category, name));
} // addJournal


// ---------------------------------------------------------------------------------------------------------------------
// Activate journal.
void
pylith::testing::TestDriver::_activateJournals(void) {
    for (size_t i = 0; i < _journals.size(); ++i) {
        JournalEnum category = _journals[i].first;
        const char* const name = _journals[i].second.c_str();
        switch (category) {
        case JOURNAL_INFO: {
            pythia::journal::info_t journal(name);
            journal.activate();
            break;
        } // INFO
        case JOURNAL_DEBUG: {
            pythia::journal::debug_t journal(name);
            journal.activate();
            break;
        } // DEBUG
        case JOURNAL_WARNING: {
            pythia::journal::warning_t journal(name);
            journal.activate();
            break;
        } // INFO
        default:
            PYLITH_JOURNAL_LOGICERROR("Unknown journal category '"<<category<<"'.");
        } // switch
    } // for
} // _activateJournal


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
// Find test matching name in test hierarchy.
const CppUnit::Test*
pylith::testing::TestDriver::_findTest(const CppUnit::Test* test,
                                       const std::string& name) {
    if (!test) { return NULL;}
    if (test->getName() == name) { return test; }
    if (!test->getChildTestCount()) { return NULL; }

    for (int i = 0; i < test->getChildTestCount(); ++i) {
        const CppUnit::Test* testCase = _findTest(test->getChildTestAt(i), name);
        if (testCase) { return testCase; }
    } // for

    return NULL;
} // _findTest


// ------------------------------------------------------------------------------------------------
int
main(int argc,
     char* argv[]) {
    pylith::testing::TestDriver driver;
    return driver.run(argc, argv);
} // main


// End of file
