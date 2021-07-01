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

#include "OutputSubfield.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/fekernels/Solution.hh" // USES Solution::passThruSubfield

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSubfield::OutputSubfield(void) :
    _dm(NULL),
    _vector(NULL),
    _fn(pylith::fekernels::Solution::passThruSubfield) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSubfield::~OutputSubfield(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSubfield::deallocate(void) {
    PetscErrorCode err;
    err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_vector);PYLITH_CHECK_ERROR(err);
} // deallocate


// ------------------------------------------------------------------------------------------------
// Create OutputSubfield from Field.
pylith::meshio::OutputSubfield*
pylith::meshio::OutputSubfield::create(const pylith::topology::Field& field,
                                       const pylith::topology::Mesh& mesh,
                                       const char* name,
                                       const int basisOrder) {
    PYLITH_METHOD_BEGIN;

    OutputSubfield* subfield = new OutputSubfield();assert(subfield);

    const pylith::topology::Field::SubfieldInfo& info = field.getSubfieldInfo(name);
    subfield->_subfieldIndex = info.index;
    subfield->_description = info.description;
    subfield->_discretization = info.fe;
    subfield->_discretization.dimension = mesh.getDimension();
    // Basis order of output should be less than or equai to the basis order of the computed field.
    subfield->_discretization.basisOrder = std::min(basisOrder, info.fe.basisOrder);

    PetscErrorCode err;
    err = DMClone(mesh.getDM(), &subfield->_dm);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_dm, name);PYLITH_CHECK_ERROR(err);

    PetscFE fe = pylith::topology::FieldOps::createFE(subfield->_discretization, subfield->_dm,
                                                      info.description.numComponents);assert(fe);
    err = PetscFESetName(fe, info.description.label.c_str());PYLITH_CHECK_ERROR(err);
    err = DMSetField(subfield->_dm, 0, NULL, (PetscObject)fe);PYLITH_CHECK_ERROR(err);
    err = DMSetFieldAvoidTensor(subfield->_dm, 0, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
    err = PetscFEDestroy(&fe);PYLITH_CHECK_ERROR(err);
    err = DMCreateDS(subfield->_dm);PYLITH_CHECK_ERROR(err);

    err = DMCreateGlobalVector(subfield->_dm, &subfield->_vector);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_vector, name);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(subfield);
}


// ------------------------------------------------------------------------------------------------
// Create OutputSubfield from Field.
pylith::meshio::OutputSubfield*
pylith::meshio::OutputSubfield::create(const pylith::topology::Field& field,
                                       const pylith::topology::Mesh& mesh,
                                       const char* name) {
    PYLITH_METHOD_BEGIN;

    OutputSubfield* subfield = new OutputSubfield();assert(subfield);

    const pylith::topology::Field::SubfieldInfo& info = field.getSubfieldInfo(name);
    subfield->_subfieldIndex = info.index;
    subfield->_description = info.description;

    PetscErrorCode err;
    err = DMClone(mesh.getDM(), &subfield->_dm);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_dm, name);PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh fieldVisitor(field, name);

    PetscSection subfieldSection = NULL;
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionClone(fieldVisitor.localSection(), &subfieldSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetChart(fieldVisitor.localSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    for (PetscInt point = pStart, offset = 0; point < pEnd; ++point) {
        const PetscInt numDof = fieldVisitor.sectionDof(point);
        err = PetscSectionSetOffset(subfieldSection, point, offset);PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetDof(subfieldSection, point, numDof);PYLITH_CHECK_ERROR(err);
        offset += numDof;
    } // for
    err = DMSetLocalSection(subfield->_dm, subfieldSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&subfieldSection);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(subfield->_dm, &subfield->_vector);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_vector, name);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(subfield);
}


// ------------------------------------------------------------------------------------------------
// Get description of subfield.
const pylith::topology::FieldBase::Description&
pylith::meshio::OutputSubfield::getDescription(void) const {
    return _description;
}


// ------------------------------------------------------------------------------------------------
// Get basis order of subfield.
int
pylith::meshio::OutputSubfield::getBasisOrder(void) const {
    return _discretization.basisOrder;
}


// ------------------------------------------------------------------------------------------------
// Get filtered PETSc global vector.
PetscVec
pylith::meshio::OutputSubfield::getVector(void) const {
    return _vector;
}


// ------------------------------------------------------------------------------------------------
// Get PETSc DM for filtered vector.
PetscDM
pylith::meshio::OutputSubfield::getDM(void) const {
    return _dm;
}


// ------------------------------------------------------------------------------------------------
// Extract subfield data from global PETSc vector with subfields.
void
pylith::meshio::OutputSubfield::project(const PetscVec& fieldVector) {
    PYLITH_METHOD_BEGIN;
    assert(fieldVector);
    assert(_vector);

    PetscErrorCode err;
    const PetscReal t = PetscReal(_subfieldIndex) + 0.01; // :KLUDGE: Easiest way to get subfield to extract into fn.
    err = DMProjectField(_dm, t, fieldVector, &_fn, INSERT_VALUES, _vector);PYLITH_CHECK_ERROR(err);
    err = VecScale(_vector, _description.scale);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Extract subfield from field.
void
pylith::meshio::OutputSubfield::extractSubfield(const pylith::topology::Field& field,
                                                const PetscInt subfieldIndex) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;
    PetscSection subfieldSection = NULL;
    PetscInt storageSize = 0;
    err = PetscSectionGetField(field.getLocalSection(), subfieldIndex, &subfieldSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetStorageSize(subfieldSection, &storageSize);PYLITH_CHECK_ERROR(err);

    PetscVec subfieldVector = this->getVector();
    PetscInt subfieldSize = 0;
    err = VecGetLocalSize(subfieldVector, &subfieldSize);PYLITH_CHECK_ERROR(err);
    assert(subfieldSize == storageSize);

    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(subfieldSection, &pStart, &pEnd);

    pylith::topology::VecVisitorMesh fieldVisitor(field);
    PetscScalar* solnArray = fieldVisitor.localArray();
    PetscScalar* subfieldArray = NULL;
    err = VecGetArray(subfieldVector, &subfieldArray);PYLITH_CHECK_ERROR(err);

    for (PetscInt point = pStart, indexVec = 0; point < pEnd; ++point) {
        const PetscInt solnOffset = fieldVisitor.sectionSubfieldOffset(subfieldIndex, point);
        const PetscInt solnDof = fieldVisitor.sectionSubfieldDof(subfieldIndex, point);

        for (PetscInt iDof = 0; iDof < solnDof; ++iDof) {
            // Dimensionalize values while extracting subfield.
            subfieldArray[indexVec++] = solnArray[solnOffset+iDof] * _description.scale;
        } // for
    } // for

    err = VecRestoreArray(subfieldVector, &subfieldArray);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // extractSubfield


// End of file
