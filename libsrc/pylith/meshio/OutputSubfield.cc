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

#include "pylith/utils/error.hh" // USES PylithCallPetsc()
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

    PylithCallPetsc(VecDestroy(&_projectVector));
    PylithCallPetsc(VecDestroy(&_projectVectorInterp));
    PylithCallPetsc(DMDestroy(&_projectDM));
    PylithCallPetsc(VecDestroy(&_outputVector));
    PylithCallPetsc(DMDestroy(&_outputDM));

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

    // Setup PETSc DM for projection
    const char* meshName = PETSC_NULLPTR;
    PylithCallPetsc(PetscObjectGetName((PetscObject) mesh.getDM(), &meshName));
    const std::string& projectName = meshName + std::string(" ") + std::string(name);
    PylithCallPetsc(DMClone(mesh.getDM(), &subfield->_projectDM));
    PylithCallPetsc(PetscObjectSetName((PetscObject)subfield->_projectDM, projectName.c_str()));
    PylithCallPetsc(DMReorderSectionSetDefault(subfield->_projectDM, DM_REORDER_DEFAULT_FALSE));
    PylithCallPetsc(DMReorderSectionSetType(subfield->_projectDM, NULL));
    PylithCallPetsc(DMPlexReorderSetDefault(subfield->_projectDM, DM_REORDER_DEFAULT_FALSE));

    // Setup PETSc FE (discretization) for projection
    PetscFE projectFE = pylith::topology::FieldOps::createFE(projectDiscretization, subfield->_projectDM,
                                                             info.description.numComponents);assert(projectFE);
    PylithCallPetsc(PetscFESetName(projectFE, info.description.label.c_str()));
    PylithCallPetsc(DMSetField(subfield->_projectDM, 0, NULL, (PetscObject)projectFE));
    PylithCallPetsc(DMSetFieldAvoidTensor(subfield->_projectDM, 0, PETSC_TRUE));
    PylithCallPetsc(PetscFEDestroy(&projectFE));
    PylithCallPetsc(DMCreateDS(subfield->_projectDM));

    if (!refineLevels) {
        subfield->_outputDM = subfield->_projectDM;
        PylithCallPetsc(PetscObjectReference((PetscObject)subfield->_outputDM));
    } else {
        delete subfield->_interpolator;subfield->_interpolator = new pylith::topology::RefineInterpolator();
        assert(subfield->_interpolator);
        subfield->_interpolator->initialize(subfield->_projectDM, refineLevels, outputBasisOrder, info.description, subfield->_discretization);
        subfield->_outputDM = subfield->_interpolator->getOutputDM();
        PylithCallPetsc(PetscObjectReference((PetscObject)subfield->_outputDM));
    } // if/else

    PylithCallPetsc(DMCreateGlobalVector(subfield->_projectDM, &subfield->_projectVector));
    PylithCallPetsc(PetscObjectSetName((PetscObject)subfield->_projectVector, name));
    if (refineLevels) {
        PylithCallPetsc(DMCreateGlobalVector(subfield->_outputDM, &subfield->_outputVector));
        PylithCallPetsc(PetscObjectSetName((PetscObject)subfield->_outputVector, name));
        PylithCallPetsc(VecDuplicate(subfield->_projectVector, &subfield->_projectVectorInterp));
    } else {
        subfield->_outputVector = subfield->_projectVector;
        PylithCallPetsc(PetscObjectReference((PetscObject)subfield->_outputVector));
    } // if/else

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

    PylithCallPetsc(DMClone(mesh.getDM(), &subfield->_projectDM));assert(subfield->_projectDM);
    PylithCallPetsc(DMReorderSectionSetDefault(subfield->_projectDM, DM_REORDER_DEFAULT_FALSE));
    PylithCallPetsc(DMReorderSectionSetType(subfield->_projectDM, NULL));
    PylithCallPetsc(PetscObjectSetName((PetscObject)subfield->_projectDM, name));

    pylith::topology::VecVisitorMesh fieldVisitor(field, name);
    fieldVisitor.selectSection(pylith::topology::VecVisitorMesh::LOCAL_SECTION);

    PetscSection subfieldSection = NULL;
    PetscInt pStart = 0, pEnd = 0;
    PylithCallPetsc(PetscSectionClone(fieldVisitor.selectedSection(), &subfieldSection));
    PylithCallPetsc(PetscSectionGetChart(fieldVisitor.selectedSection(), &pStart, &pEnd));
    for (PetscInt point = pStart, offset = 0; point < pEnd; ++point) {
        const PetscInt numDof = fieldVisitor.sectionDof(point);
        PylithCallPetsc(PetscSectionSetOffset(subfieldSection, point, offset));
        PylithCallPetsc(PetscSectionSetDof(subfieldSection, point, numDof));
        offset += numDof;
    } // for
    PylithCallPetsc(DMSetLocalSection(subfield->_projectDM, subfieldSection));
    PylithCallPetsc(PetscSectionDestroy(&subfieldSection));
    PylithCallPetsc(DMCreateGlobalVector(subfield->_projectDM, &subfield->_projectVector));
    PylithCallPetsc(PetscObjectSetName((PetscObject)subfield->_projectVector, name));

    subfield->_outputDM = subfield->_projectDM;
    PylithCallPetsc(PetscObjectReference((PetscObject) subfield->_outputDM));
    subfield->_outputVector = subfield->_projectVector;
    PylithCallPetsc(PetscObjectReference((PetscObject) subfield->_outputVector));

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
    PylithCallPetsc(DMGetLabel(_projectDM, name, &_label));
    PylithCallPetsc(DMPlexLabelComplete(_projectDM, _label));

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

    const PetscReal t = PetscReal(_subfieldIndex) + 0.01; // :KLUDGE: Easiest way to get subfield to extract into fn

    PylithCallPetsc(DMProjectField(_projectDM, t, fieldVector, &_fn, INSERT_VALUES, _projectVector));
    if (_interpolator) {
        pylith::topology::FieldOps::transformVector(&_projectVectorInterp, _interpolator->getInputDM(), _projectVector, _projectDM);
        _interpolator->interpolate(&_outputVector, _projectVectorInterp);
    } // if
    PylithCallPetsc(VecScale(_outputVector, _description.scale));

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

    const PetscReal t = PetscReal(_subfieldIndex) + 0.01; // :KLUDGE: Easiest way to get subfield to extract into fn

    PylithCallPetsc(DMProjectFieldLabel(_projectDM, t, _label, 1, &_labelValue, PETSC_DETERMINE, NULL, fieldVector, &_fn, INSERT_VALUES, _projectVector));
    if (_interpolator) {
        pylith::topology::FieldOps::transformVector(&_projectVectorInterp, _interpolator->getInputDM(), _projectVector, _projectDM);
        _interpolator->interpolate(&_outputVector, _projectVectorInterp);
    } // if
    PylithCallPetsc(VecScale(_outputVector, _description.scale));

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

    PetscSection subfieldSection = NULL;
    PetscInt storageSize = 0;
    PylithCallPetsc(PetscSectionGetField(field.getGlobalSection(), subfieldIndex, &subfieldSection));
    PylithCallPetsc(PetscSectionGetStorageSize(subfieldSection, &storageSize));

    PetscVec subfieldVector = this->getOutputVector();
    PetscInt subfieldSize = 0;
    PylithCallPetsc(VecGetLocalSize(subfieldVector, &subfieldSize));
    assert(subfieldSize == storageSize);

    assert(field.getOutputVector());
    PetscIS subfieldIS = PETSC_NULLPTR;
    PylithCallPetsc(DMCreateSubDM(field.getDM(), 1, &subfieldIndex, &subfieldIS, PETSC_NULLPTR));
    PylithCallPetsc(VecISCopy(field.getOutputVector(), subfieldIS, SCATTER_REVERSE, subfieldVector));
    PylithCallPetsc(VecScale(subfieldVector, _description.scale));
    PylithCallPetsc(ISDestroy(&subfieldIS));

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::extractSubfield);
    PYLITH_METHOD_END;
} // extractSubfield


// End of file
