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

#include "StubMethodTracker.hh" // Implementation of class methods

pylith::testing::StubMethodTracker::map_type pylith::testing::StubMethodTracker::_methodCount;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::testing::StubMethodTracker::StubMethodTracker(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add to count for method.
pylith::testing::StubMethodTracker::StubMethodTracker(const char* methodName) {
    methodCalled(methodName);
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::testing::StubMethodTracker::~StubMethodTracker(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add to count for method.
void
pylith::testing::StubMethodTracker::methodCalled(const char* methodName) {
    map_type::iterator iter = _methodCount.find(std::string(methodName));
    if (iter != _methodCount.end()) {
        iter->second++;
    } else {
        _methodCount[std::string(methodName)] = 1;
    } // if/else
} // methodCalled


// ---------------------------------------------------------------------------------------------------------------------
// Reset method counts.
void
pylith::testing::StubMethodTracker::clear(void) {
    _methodCount.clear();
} // clear


// ---------------------------------------------------------------------------------------------------------------------
// How many times was method called?
size_t
pylith::testing::StubMethodTracker::getMethodCount(const char* methodName) {
    size_t count = 0;

    const map_type::iterator iter = _methodCount.find(std::string(methodName));
    if (iter != _methodCount.end()) {
        count = iter->second;
    } // if

    return count;
} // getMethodCount


// End of file
