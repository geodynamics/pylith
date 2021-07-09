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

/**
 * @file libsrc/problems/TestDriver.hh
 *
 * @brief Object for running C++ tests.
 */

#if !defined(pylith_testing_testdriver_hh)
#define pylith_testing_testdriver_hh

#include "pylith/testing/testingfwd.hh" // forward declaration

#include <cppunit/Test.h> /// USES CppUnit::Test

#include <vector> // USES std::vector
#include <utility> // USES std::pair
#include <string> // USES std::string

class pylith::testing::TestDriver {
public:

    enum JournalEnum {
        JOURNAL_INFO=0,
        JOURNAL_DEBUG=1,
        JOURNAL_WARNING=2,
    };

    typedef std::vector< std::pair<JournalEnum,std::string> > journals_t;

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestDriver(void);

    /// Destructor
    ~TestDriver(void);

    /** Run test application.
     *
     * Arguments:
     *   --help
     *   --list
     *   --quiet
     *   --tests=TEST_0,...,TEST_N
     *   --petsc OPTION=VALUE
     *   --journal.info=NAME
     *   --journal.debug=NAME
     *   --journal.warning=NAME
     *
     * @param argc[in] Number of arguments passed.
     * @param argv[in] Array of input arguments.
     *
     * @returns 1 if errors were detected, 0 otherwise.
     */
    int run(int argc,
            char* argv[]);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
    /** Parse command line arguments.
     *
     * @param[in] argc Number of arguments passed.
     * @param[in] argv Array of input arguments.
     */
    void _parseArgs(int argc,
                    char* argv[]);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::vector<std::string> _tests;
    std::vector<std::string> _petscOptions;
    journals_t _journals;
    bool _showHelp;
    bool _listTests;
    bool _mallocDump;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    TestDriver(const TestDriver&); ///< Not implemented
    const TestDriver& operator=(const TestDriver&); ///< Not implemented

}; // TestDriver

#endif // pylith_testing_testdriver_hh

// End of file
