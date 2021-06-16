// -*- C++ -*-
//
// ---------------------------------------------------------------------------------------------------------------------
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
// ---------------------------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "Source.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::Source::Source(void) :
    _sourceId(0),
    _descriptiveLabel("") {
    //
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::Source::~Source(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::sources::Source::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::problems::Physics::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set value of label source-id used to identify source cells.
void
pylith::sources::Source::setSourceId(const int value) {
    PYLITH_COMPONENT_DEBUG("setsourceId(value="<<value<<")");

    _sourceId = value;
} // setSourceId


// ---------------------------------------------------------------------------------------------------------------------
// Get value of label source-id used to identify source cells.
int
pylith::sources::Source::getSourceId(void) const {
    return _sourceId;
} // getSourceId


// ---------------------------------------------------------------------------------------------------------------------
// Set descriptive label of source.
void
pylith::sources::Source::setDescriptiveLabel(const char* value) {
    PYLITH_COMPONENT_DEBUG("setDescriptiveLabel(value="<<value<<")");

    _descriptiveLabel = value;
} // setDescriptiveLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label of source.
const char*
pylith::sources::Source::getDescriptiveLabel(void) const {
    return _descriptiveLabel.c_str();
} // getDescriptiveLabel


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::sources::Source::createConstraint(const pylith::topology::Field& solution) {
    return NULL;
} // createConstraint


// ------------------------------------------------------------------------------------------------
// Set point names and coordinates of points .
void
pylith::sources::Source::setPoints(const PylithReal* pointCoords,
                                   const int numPoints,
                                   const int spaceDim,
                                   const char* const* pointNames,
                                   const int numPointNames) {
    PYLITH_METHOD_BEGIN;

    assert(pointCoords && pointNames);
    assert(numPoints == numPointNames);

    // Copy point coordinates.
    const PylithInt size = numPoints * spaceDim;
    _pointCoords.resize(size);
    for (PylithInt i = 0; i < size; ++i) {
        _pointCoords[i] = pointCoords[i];
    } // for

    // Copy point names.
    _pointNames.resize(numPointNames);
    for (PylithInt i = 0; i < numPointNames; ++i) {
        _pointNames[i] = pointNames[i];
    } // for

    PYLITH_METHOD_END;
} // setPoints


// End of file
