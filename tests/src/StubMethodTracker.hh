// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

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

// End of file
