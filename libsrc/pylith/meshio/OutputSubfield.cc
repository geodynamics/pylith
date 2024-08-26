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

#include "pylith/meshio/OutputSubfield.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/RefineInterpolator.hh" // USES RefineInterpolator
#include "pylith/fekernels/Solution.hh" // USES Solution::passThruSubfield

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "petscdm.h" // USES DMReorderSectionSetDefault()

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace meshio {
        class _OutputSubfield {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt create;
                static PylithInt createBasisOrder;
                static PylithInt setLabel;
                static PylithInt project;
                static PylithInt projectWithLabel;
                static PylithInt extractSubfield;
            };

        }; // _OutputSubfield
    } // meshio
} // pylith

pylith::utils::EventLogger pylith::meshio::_OutputSubfield::Events::logger;
PylithInt pylith::meshio::_OutputSubfield::Events::create;
PylithInt pylith::meshio::_OutputSubfield::Events::createBasisOrder;
PylithInt pylith::meshio::_OutputSubfield::Events::setLabel;
PylithInt pylith::meshio::_OutputSubfield::Events::project;
PylithInt pylith::meshio::_OutputSubfield::Events::projectWithLabel;
PylithInt pylith::meshio::_OutputSubfield::Events::extractSubfield;

// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_OutputSubfield::Events::init(void) {
    logger.setClassName("OutputSubfield");
    logger.initialize();
    create = logger.registerEvent("PL:OutputSubfield:create");
    createBasisOrder = logger.registerEvent("PL:OutputSubfield:createBasisOrder");
    setLabel = logger.registerEvent("PL:OutputSubfield:setLabel");
    project = logger.registerEvent("PL:OutputSubfield:project");
    projectWithLabel = logger.registerEvent("PL:OutputSubfield:projectWithLabel");
    extractSubfield = logger.registerEvent("PL:OutputSubfield:extractSubfield");
}


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSubfield::OutputSubfield(void) :
    _subfieldIndex(-1),
    _projectDM(PETSC_NULLPTR),
    _projectVector(PETSC_NULLPTR),
    _projectVectorInterp(PETSC_NULLPTR),
    _fn(pylith::fekernels::Solution::passThruSubfield),
    _outputDM(PETSC_NULLPTR),
    _outputVector(PETSC_NULLPTR),
    _interpolator(NULL),
    _label(NULL),
    _labelValue(0) {
    _OutputSubfield::Events::init();
}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSubfield::~OutputSubfield(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSubfield::deallocate(void) {
    delete _interpolator;_interpolator = NULL;

    PetscErrorCode err;
    err = VecDestroy(&_projectVector);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_projectVectorInterp);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&_projectDM);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_outputVector);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&_outputDM);PYLITH_CHECK_ERROR(err);

    _label = NULL; // Destroyed by DMDestroy()
} // deallocate


// ------------------------------------------------------------------------------------------------
// Create OutputSubfield from Field.
pylith::meshio::OutputSubfield*
pylith::meshio::OutputSubfield::create(const pylith::topology::Field& field,
                                       const pylith::topology::Mesh& mesh,
                                       const char* name,
                                       const int basisOrder,
                                       const int refineLevels) {
    PYLITH_METHOD_BEGIN;
    // _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::create);

    OutputSubfield* subfield = new OutputSubfield();assert(subfield);

    const pylith::topology::Field::SubfieldInfo& info = field.getSubfieldInfo(name);
    subfield->_subfieldIndex = info.index;
    subfield->_description = info.description;
    const int outputBasisOrder = std::min(basisOrder, info.fe.basisOrder);

    // Discretization for projection
    pylith::topology::FieldBase::Discretization projectDiscretization = info.fe;
    projectDiscretization.dimension = mesh.getDimension();
    projectDiscretization.basisOrder = refineLevels ? info.fe.basisOrder : outputBasisOrder;

    // Discretization for output
    subfield->_discretization = projectDiscretization;
    subfield->_discretization.basisOrder = outputBasisOrder;

    PetscErrorCode err = PETSC_SUCCESS;

    // Setup PETSc DM for projection
    err = DMClone(mesh.getDM(), &subfield->_projectDM);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_projectDM, name);PYLITH_CHECK_ERROR(err);
    err = DMReorderSectionSetDefault(subfield->_projectDM, DM_REORDER_DEFAULT_FALSE);PYLITH_CHECK_ERROR(err);
    err = DMReorderSectionSetType(subfield->_projectDM, NULL);PYLITH_CHECK_ERROR(err);
    err = DMPlexReorderSetDefault(subfield->_projectDM, DM_REORDER_DEFAULT_FALSE);

    // Setup PETSc FE (discretization) for projection
    PetscFE projectFE = pylith::topology::FieldOps::createFE(projectDiscretization, subfield->_projectDM,
                                                             info.description.numComponents);assert(projectFE);
    err = PetscFESetName(projectFE, info.description.label.c_str());PYLITH_CHECK_ERROR(err);
    err = DMSetField(subfield->_projectDM, 0, NULL, (PetscObject)projectFE);PYLITH_CHECK_ERROR(err);
    err = DMSetFieldAvoidTensor(subfield->_projectDM, 0, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
    err = PetscFEDestroy(&projectFE);PYLITH_CHECK_ERROR(err);
    err = DMCreateDS(subfield->_projectDM);PYLITH_CHECK_ERROR(err);

    if (!refineLevels) {
        subfield->_outputDM = subfield->_projectDM;
        err = PetscObjectReference((PetscObject)subfield->_outputDM);PYLITH_CHECK_ERROR(err);
    } else {
        delete subfield->_interpolator;subfield->_interpolator = new pylith::topology::RefineInterpolator();
        assert(subfield->_interpolator);
        subfield->_interpolator->initialize(subfield->_projectDM, refineLevels, outputBasisOrder, info.description, subfield->_discretization);
        subfield->_outputDM = subfield->_interpolator->getOutputDM();
        err = PetscObjectReference((PetscObject)subfield->_outputDM);PYLITH_CHECK_ERROR(err);
    } // if/else

    err = DMCreateGlobalVector(subfield->_projectDM, &subfield->_projectVector);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_projectVector, name);PYLITH_CHECK_ERROR(err);
    if (refineLevels) {
        err = DMCreateGlobalVector(subfield->_outputDM, &subfield->_outputVector);PYLITH_CHECK_ERROR(err);
        err = PetscObjectSetName((PetscObject)subfield->_outputVector, name);PYLITH_CHECK_ERROR(err);
        err = VecDuplicate(subfield->_projectVector, &subfield->_projectVectorInterp);PYLITH_CHECK_ERROR(err);
    } else {
        subfield->_outputVector = subfield->_projectVector;
        err = PetscObjectReference((PetscObject)subfield->_outputVector);PYLITH_CHECK_ERROR(err);
    }

    // _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::create);
    PYLITH_METHOD_RETURN(subfield);
}


// ------------------------------------------------------------------------------------------------
// Create OutputSubfield from Field.
pylith::meshio::OutputSubfield*
pylith::meshio::OutputSubfield::create(const pylith::topology::Field& field,
                                       const pylith::topology::Mesh& mesh,
                                       const char* name) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::createBasisOrder);

    OutputSubfield* subfield = new OutputSubfield();assert(subfield);

    const pylith::topology::Field::SubfieldInfo& info = field.getSubfieldInfo(name);
    subfield->_subfieldIndex = info.index;
    subfield->_description = info.description;

    PetscErrorCode err = PETSC_SUCCESS;
    err = DMClone(mesh.getDM(), &subfield->_projectDM);PYLITH_CHECK_ERROR(err);
    err = DMReorderSectionSetDefault(subfield->_projectDM, DM_REORDER_DEFAULT_FALSE);PYLITH_CHECK_ERROR(err);
    err = DMReorderSectionSetType(subfield->_projectDM, NULL);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_projectDM, name);PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh fieldVisitor(field, name);

    PetscSection subfieldSection = NULL;
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionClone(fieldVisitor.selectedSection(), &subfieldSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetChart(fieldVisitor.selectedSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    for (PetscInt point = pStart, offset = 0; point < pEnd; ++point) {
        const PetscInt numDof = fieldVisitor.sectionDof(point);
        err = PetscSectionSetOffset(subfieldSection, point, offset);PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetDof(subfieldSection, point, numDof);PYLITH_CHECK_ERROR(err);
        offset += numDof;
    } // for
    err = DMSetLocalSection(subfield->_projectDM, subfieldSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&subfieldSection);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(subfield->_projectDM, &subfield->_projectVector);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_projectVector, name);PYLITH_CHECK_ERROR(err);

    subfield->_outputDM = subfield->_projectDM;
    err = PetscObjectReference((PetscObject) subfield->_outputDM);PYLITH_CHECK_ERROR(err);
    subfield->_outputVector = subfield->_projectVector;
    err = PetscObjectReference((PetscObject) subfield->_outputVector);PYLITH_CHECK_ERROR(err);

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::createBasisOrder);
    PYLITH_METHOD_RETURN(subfield);
}


// ------------------------------------------------------------------------------------------------
// Set label name and value.
void
pylith::meshio::OutputSubfield::setLabel(const char* name,
                                         const int value) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::setLabel);

    if (_label) {
        PYLITH_METHOD_END;
    } // if
    PetscErrorCode err;
    err = DMGetLabel(_projectDM, name, &_label);PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelComplete(_projectDM, _label);

    _labelValue = value;

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::setLabel);
    PYLITH_METHOD_END;
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
pylith::meshio::OutputSubfield::getOutputVector(void) const {
    return _outputVector;
}


// ------------------------------------------------------------------------------------------------
// Get PETSc DM for filtered vector.
PetscDM
pylith::meshio::OutputSubfield::getOutputDM(void) const {
    return _outputDM;
}


// ------------------------------------------------------------------------------------------------
// Extract subfield data from global PETSc vector with subfields.
void
pylith::meshio::OutputSubfield::project(const PetscVec& fieldVector) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::project);
    assert(fieldVector);
    assert(_projectVector);
    assert(_outputVector);

    PetscErrorCode err;
    const PetscReal t = PetscReal(_subfieldIndex) + 0.01; // :KLUDGE: Easiest way to get subfield to extract into fn

    err = DMProjectField(_projectDM, t, fieldVector, &_fn, INSERT_VALUES, _projectVector);PYLITH_CHECK_ERROR(err);
    if (_interpolator) {
        pylith::topology::FieldOps::transformVector(&_projectVectorInterp, _interpolator->getInputDM(), _projectVector, _projectDM);
        _interpolator->interpolate(&_outputVector, _projectVectorInterp);
    } // if
    err = VecScale(_outputVector, _description.scale);PYLITH_CHECK_ERROR(err);

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::project);
    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Extract subfield data from global PETSc vector with subfields.
void
pylith::meshio::OutputSubfield::projectWithLabel(const PetscVec& fieldVector) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::projectWithLabel);
    assert(fieldVector);
    assert(_projectVector);
    assert(_outputVector);
    assert(_label);

    PetscErrorCode err;
    const PetscReal t = PetscReal(_subfieldIndex) + 0.01; // :KLUDGE: Easiest way to get subfield to extract into fn

    err = DMProjectFieldLabel(_projectDM, t, _label, 1, &_labelValue, PETSC_DETERMINE, NULL, fieldVector, &_fn, INSERT_VALUES, _projectVector);PYLITH_CHECK_ERROR(err);
    if (_interpolator) {
        pylith::topology::FieldOps::transformVector(&_projectVectorInterp, _interpolator->getInputDM(), _projectVector, _projectDM);
        _interpolator->interpolate(&_outputVector, _projectVectorInterp);
    } // if
    err = VecScale(_outputVector, _description.scale);PYLITH_CHECK_ERROR(err);

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::projectWithLabel);
    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Extract subfield from field.
void
pylith::meshio::OutputSubfield::extractSubfield(const pylith::topology::Field& field,
                                                const PetscInt subfieldIndex) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::extractSubfield);

    PetscErrorCode err;
    PetscSection subfieldSection = NULL;
    PetscInt storageSize = 0;
    err = PetscSectionGetField(field.getLocalSection(), subfieldIndex, &subfieldSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetStorageSize(subfieldSection, &storageSize);PYLITH_CHECK_ERROR(err);

    PetscVec subfieldVector = this->getOutputVector();
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

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::extractSubfield);
    PYLITH_METHOD_END;
} // extractSubfield


// End of file
