// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#include <portinfo>

#include "pylith/testing/testingfwd.hh"
#include "pylith/topology/FieldOps.hh" // USES FieldOps::deallocate()
#include "pylith/utils/journals.hh" // USES journals

#include "catch2/catch_session.hpp"

#include "petsc.h"
#include <Python.h>

#include <getopt.h> // USES getopt_long()
#include <sstream> // USES std::ostringstream, std::istringstream
#include <vector> // USES std::vector
#include <utility> // USES std::pair
#include <string> // USES std::string
#include <iostream> // USES std::cout

// ------------------------------------------------------------------------------------------------
class pylith::testing::TestDriver {
public:

    enum JournalEnum {
        JOURNAL_INFO=0,
        JOURNAL_DEBUG=1,
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

    /** Initialize PETSc.
     *
     * @param[in] programName Name of executable program.
     * @param[in] petscArgs String with comma separated list of PETSc arguments.
     * @param[in] mallocDump Set malloc debug dump.
     */
    int _initializePetsc(char* programName,
                         const std::string& petscArgs,
                         const bool mallocDump);

    /** Activate journals.
     *
     * @param[in] category Journal category.
     * @param[in] journalArgs String with comma separted list of journal names.
     */
    void _activateJournals(JournalEnum category,
                           const std::string& journalArgs);

    /** Extract values from string delimited by a character.
     *
     * @param[in] values String with delimiter separated list of values.
     * @param[in] delimiter Character delimiting values.
     */
    std::vector<std::string> _parseList(const std::string& values,
                                        const char delimiter=',');

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    TestDriver(const TestDriver&); ///< Not implemented
    const TestDriver& operator=(const TestDriver&); ///< Not implemented

};

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::testing::TestDriver::TestDriver(void) { }


// ---------------------------------------------------------------------------------------------------------------------
// Run info application.
int
pylith::testing::TestDriver::run(int argc,
                                 char* argv[]) {
    Catch::Session session;

    bool quiet = false;
    std::string petscArgs;
    std::string infoJournalArgs;
    std::string debugJournalArgs;
    auto cli = session.cli()
               | Catch::Clara::Opt(quiet, "quiet")["--petsc-quiet"]("Turn off some PETSc debugging flags (e.g., memory leak detection)")
               | Catch::Clara::Opt(petscArgs, "petsc")["--petsc"]("Comma separated list of PETSc arguments")
               | Catch::Clara::Opt(infoJournalArgs, "journal.info")["--journal.info"]("Comma separated list of info journal names to activate")
               | Catch::Clara::Opt(debugJournalArgs, "journal.debug")["--journal.debug"]("Comma separated list of debug journal names to activate");
    session.cli(cli);
    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode) {
        return returnCode;
    } // if

    // Initialize PETSc
#if defined(MALLOC_DEBUG_OFF)
    const bool mallocDump = false;
#else
    const bool mallocDump = !quiet;
#endif
    int err = _initializePetsc(argv[0], petscArgs, mallocDump);CHKERRQ(err);

    // Initialize Python (needed for journals).
    Py_Initialize();

    _activateJournals(JOURNAL_INFO, infoJournalArgs);
    _activateJournals(JOURNAL_DEBUG, debugJournalArgs);
    int result = session.run();

    // Finalize Python
    Py_Finalize();

    // Finalize PETSc
    pylith::topology::FieldOps::deallocate();
    err = PetscFinalize();CHKERRQ(err);

    if (!mallocDump) {
        std::cout << "WARNING PETSc option -malloc dump is OFF\n" << std::endl;
    } // if

    return result;
} // run


// ---------------------------------------------------------------------------------------------------------------------
// Initialize PETSc.
int
pylith::testing::TestDriver::_initializePetsc(char* programName,
                                              const std::string& petscArgs,
                                              const bool mallocDump) {
    int argc = 1;
    char** argv = new char*[argc+1];
    argv[0] = programName;
    argv[argc] = nullptr; // C standard is argv[argc] == nullptr.
    PetscErrorCode err = PetscInitialize(&argc, &argv, nullptr, nullptr);CHKERRQ(err);
    delete[] argv;argv = nullptr;

    if (mallocDump) {
        err = PetscOptionsSetValue(nullptr, "-malloc_dump", "");CHKERRQ(err);
    } // if

    const std::vector<std::string>& petscOptions = _parseList(petscArgs);
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


// ------------------------------------------------------------------------------------------------
// Activate journal.
void
pylith::testing::TestDriver::_activateJournals(JournalEnum category,
                                               const std::string& journalArgs) {
    const std::vector<std::string>& names = _parseList(journalArgs);
    for (size_t i = 0; i < names.size(); ++i) {
        const std::string& name = names[i];
        switch (category) {
        case JOURNAL_INFO: {
            pythia::journal::info_t(name).activate();
            break;
        } // INFO
        case JOURNAL_DEBUG: {
            pythia::journal::debug_t(name).activate();
            break;
        } // DEBUG
        default:
            PYLITH_JOURNAL_LOGICERROR("Unknown journal category '"<<category<<"'.");
        } // switch
    } // for
} // _activateJournal


// ------------------------------------------------------------------------------------------------
// Parse string delimited by a character.
std::vector<std::string>
pylith::testing::TestDriver::_parseList(const std::string& values,
                                        const char delimiter) {
    std::istringstream stream(values);
    std::vector<std::string> tokens;

    while (stream.good()) {
        std::string token;
        std::getline(stream, token, delimiter);
        if (!token.empty()) {
            tokens.push_back(token);
        } // if
    } // while

    return tokens;
} // parseList


// ------------------------------------------------------------------------------------------------
int
main(int argc,
     char* argv[]) {
    return pylith::testing::TestDriver().run(argc, argv);
} // main


// End of file
