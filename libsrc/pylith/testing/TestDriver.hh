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
#include <string> // USES std::string

class pylith::testing::TestDriver {
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
     *
     * @param argc[in] Number of arguments passed.
     * @param argv[in] Array of input arguments.
     *
     * @returns 1 if errors were detected, 0 otherwise.
     */
    int run(int argc,
            char* argv[]);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /** Parse command line arguments.
     *
     * @param argc[in] Number of arguments passed.
     * @param argv[in] Array of input arguments.
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

    /** Find test matching name in test hierarchy.
     *
     * @param[in] test Test hierarchy.
     * @param[in] name Name of test to find.
     * @returns Test matching name or NULL if not found.
     */
    const CppUnit::Test* _findTest(const CppUnit::Test* const test,
                                   const std::string& name);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::vector<std::string> _tests;
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
