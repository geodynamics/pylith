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
 * @file libsrc/testing/StubMethodTracker.hh
 *
 * @brief C++ implementation to verify which Stub methods are executed in testing.
 */

#if !defined(pylith_testing_stubmethodtracker_hh)
#define pylith_testing_stubmethodtracker_hh

#include "pylith/testing/testingfwd.hh" // forward declarations

#include <string> // USES std::string
#include <map> // USES std::map

class pylith::testing::StubMethodTracker {
    friend class TestStubMethodTracker; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    StubMethodTracker(void);

    /** Constructor with method name.
     *
     * @param[in] methodName Full namespace method name.
     */
    StubMethodTracker(const char* name);

    /// Destructor
    ~StubMethodTracker(void);

    /** Add to count for method.
     *
     * @param[in] methodName Full namespace method name.
     */
    void methodCalled(const char* methodName);

    /// Reset method counts.
    void clear(void);

    /** How many times was method called?
     *
     * @param[in] methodName Full namespace method name.
     * @returns Number of times method was called.
     */
    size_t getMethodCount(const char* methodName);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    StubMethodTracker(const StubMethodTracker&); ///< Not implemented.
    const StubMethodTracker& operator=(const StubMethodTracker&); ///< Not implemented

    typedef std::map<std::string, int> map_type;
    static map_type _methodCount; ///< Number of times method was called.

}; // StubMethodTracker

#endif // pylith_testing_stubmethodtrack_hh

// End of file
