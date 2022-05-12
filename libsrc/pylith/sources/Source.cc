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
    _description(""),
    _labelName(""),
    _labelValue(1) {
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

// ------------------------------------------------------------------------------------------------
// Set descriptive label of source.
void
pylith::sources::Source::setDescription(const char* value) {
    PYLITH_COMPONENT_DEBUG("setDescription(value="<<value<<")");

    _description = value;
} // setDescription


// ------------------------------------------------------------------------------------------------
// Get label of source.
const char*
pylith::sources::Source::getDescription(void) const {
    return _description.c_str();
} // getDescription


// ------------------------------------------------------------------------------------------------
// Set name of label marking material.
void
pylith::sources::Source::setLabelName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for material label.");
    } // if

    _labelName = value;
} // setLabelName


// ------------------------------------------------------------------------------------------------
// Get name of label marking material.
const char*
pylith::sources::Source::getLabelName(void) const {
    return _labelName.c_str();
} // getLabelName


// ------------------------------------------------------------------------------------------------
// Set value of label marking material.
void
pylith::sources::Source::setLabelValue(const int value) {
    _labelValue = value;
} // setLabelValue


// ------------------------------------------------------------------------------------------------
// Get value of label marking material.
int
pylith::sources::Source::getLabelValue(void) const {
    return _labelValue;
} // getLabelValue


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<pylith::feassemble::Constraint*>
pylith::sources::Source::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getLabel()<<") empty method");
    std::vector<pylith::feassemble::Constraint*> constraintArray;

    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


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
