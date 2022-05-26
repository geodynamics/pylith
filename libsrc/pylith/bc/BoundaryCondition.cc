// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "BoundaryCondition.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cstring> // USES strlen()
#include <stdexcept> // USES std::runtime_error()

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::BoundaryCondition::BoundaryCondition(void) :
    _subfieldName(""),
    _labelName(""),
    _labelValue(1) {
    _refDir1[0] = 0.0;
    _refDir1[1] = 0.0;
    _refDir1[2] = 1.0;

    _refDir2[0] = 0.0;
    _refDir2[1] = 1.0;
    _refDir2[2] = 0.0;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::BoundaryCondition::~BoundaryCondition(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::BoundaryCondition::deallocate(void) {
    Physics::deallocate();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of solution subfield associated with boundary condition.
void
pylith::bc::BoundaryCondition::setSubfieldName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setSubfieldName(value="<<value<<")");

    if (!value || (0 == strlen(value))) {
        std::ostringstream msg;
        msg << "Empty string given for name of solution subfield for boundary condition '" << _labelName
            <<"'.";
        throw std::runtime_error(msg.str());
    } // if
    _subfieldName = value;
} // setSubfieldName


// ------------------------------------------------------------------------------------------------
// Get name of solution subfield associated with boundary condition.
const char*
pylith::bc::BoundaryCondition::getSubfieldName(void) const {
    return _subfieldName.c_str();
} // getSubfieldName


// ------------------------------------------------------------------------------------------------
// Set name of label marking boundary associated with boundary condition surface.
void
pylith::bc::BoundaryCondition::setLabelName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition label.");
    } // if

    _labelName = value;
} // setLabelName


// ------------------------------------------------------------------------------------------------
// Get name of label marking boundary associated with boundary condition surface.
const char*
pylith::bc::BoundaryCondition::getLabelName(void) const {
    return _labelName.c_str();
} // getLabelName


// ------------------------------------------------------------------------------------------------
// Set value of label marking boundary associated with boundary condition surface.
void
pylith::bc::BoundaryCondition::setLabelValue(const int value) {
    _labelValue = value;
} // setLabelValue


// ------------------------------------------------------------------------------------------------
// Get value of label marking boundary associated with boundary condition surface.
int
pylith::bc::BoundaryCondition::getLabelValue(void) const {
    return _labelValue;
} // getLabelValue


// ------------------------------------------------------------------------------------------------
// Set first choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::bc::BoundaryCondition::setRefDir1(const PylithReal vec[3]) {
    PYLITH_COMPONENT_DEBUG("setRefDir1(vec="<<vec[0]<<","<<vec[1]<<","<<vec[2]<<")");

    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 1 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for boundary condition '" << _labelName << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    } // for
} // setRefDir1


// ------------------------------------------------------------------------------------------------
// Set second choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::bc::BoundaryCondition::setRefDir2(const PylithReal vec[3]) {
    PYLITH_COMPONENT_DEBUG("setRefDir2(vec="<<vec[0]<<","<<vec[1]<<","<<vec[2]<<")");

    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 2 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for boundary condition '" << _labelName << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir2[i] = vec[i] / mag;
    } // for
} // setRefDir2


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::BoundaryCondition::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    if (!solution.hasSubfield(_subfieldName.c_str())) {
        std::ostringstream msg;
        msg << "Cannot apply boundary condition to field '"<< _subfieldName
            << "'; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    const PetscDM dmSoln = solution.getDM();
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err = DMHasLabel(dmSoln, _labelName.c_str(), &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << _labelName << "' in boundary condition '"
            << PyreComponent::getIdentifier() <<"'.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// End of file
