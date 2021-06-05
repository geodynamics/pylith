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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Field.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh
#include "FieldOps.hh" // USES FieldOps
#include "VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/topology/MeshOps.hh" // USES isCohesiveCell()

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include <iostream>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::topology::Field::Field(const Mesh& mesh) :
    _mesh(mesh),
    _dm(NULL),
    _localVec(NULL),
    _globalVec(NULL),
    _outputVec(NULL) {
    PYLITH_METHOD_BEGIN;

    GenericComponent::setName("field");
    _label = "unknown";

    if (mesh.dmMesh()) {
        PetscDM dm = mesh.dmMesh();assert(dm);
        PetscVec coordVec = NULL;
        PetscErrorCode err;

        err = DMClone(dm, &_dm);PYLITH_CHECK_ERROR(err);
        err = DMGetCoordinatesLocal(dm, &coordVec);PYLITH_CHECK_ERROR(err);
        if (coordVec) {
            PetscDM coordDM = NULL, newCoordDM = NULL;
            PetscSection coordSection = NULL, newCoordSection = NULL;

            err = DMGetCoordinateDM(dm, &coordDM);PYLITH_CHECK_ERROR(err);
            err = DMGetCoordinateDM(_dm, &newCoordDM);PYLITH_CHECK_ERROR(err);
            err = DMGetSection(coordDM, &coordSection);PYLITH_CHECK_ERROR(err);
            err = PetscSectionClone(coordSection, &newCoordSection);PYLITH_CHECK_ERROR(err);
            err = DMSetSection(newCoordDM, newCoordSection);PYLITH_CHECK_ERROR(err);
            err = PetscSectionDestroy(&newCoordSection);PYLITH_CHECK_ERROR(err);
            err = DMSetCoordinatesLocal(_dm, coordVec);PYLITH_CHECK_ERROR(err);
        } // if

    } // if

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Constructor with field to use for layout.
pylith::topology::Field::Field(const Field& src) :
    _mesh(src._mesh),
    _dm(NULL),
    _localVec(NULL),
    _globalVec(NULL),
    _outputVec(NULL) {
    PYLITH_METHOD_BEGIN;

    GenericComponent::setName("field");
    _label = "unknown";

    PetscErrorCode err;
    const char* name = NULL;

    if (src.dmMesh()) {
        PetscDM dm = src.dmMesh();assert(dm);
        PetscVec coordVec = NULL;
        PetscErrorCode err;

        err = DMClone(dm, &_dm);PYLITH_CHECK_ERROR(err);
        err = DMGetCoordinatesLocal(dm, &coordVec);PYLITH_CHECK_ERROR(err);
        if (coordVec) {
            PetscDM coordDM = NULL, newCoordDM = NULL;
            PetscSection coordSection = NULL, newCoordSection = NULL;

            err = DMGetCoordinateDM(dm, &coordDM);PYLITH_CHECK_ERROR(err);
            err = DMGetCoordinateDM(_dm, &newCoordDM);PYLITH_CHECK_ERROR(err);
            err = DMGetSection(coordDM, &coordSection);PYLITH_CHECK_ERROR(err);
            err = PetscSectionClone(coordSection, &newCoordSection);PYLITH_CHECK_ERROR(err);
            err = DMSetSection(newCoordDM, newCoordSection);PYLITH_CHECK_ERROR(err);
            err = PetscSectionDestroy(&newCoordSection);PYLITH_CHECK_ERROR(err);
            err = DMSetCoordinatesLocal(_dm, coordVec);PYLITH_CHECK_ERROR(err);
        } // if
        err = DMCopyDisc(dm, _dm);PYLITH_CHECK_ERROR(err);
    } // if

    assert(_dm);
    PetscSection srcSection = src.localSection();
    PetscSection newSection = NULL;
    err = PetscSectionClone(srcSection, &newSection);PYLITH_CHECK_ERROR(err);
    err = DMSetSection(_dm, newSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&newSection);PYLITH_CHECK_ERROR(err);

    assert(!_localVec);
    err = DMCreateLocalVector(_dm, &_localVec);PYLITH_CHECK_ERROR(err);
    err = PetscObjectGetName((PetscObject)src._dm, &name);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)_dm,  name);PYLITH_CHECK_ERROR(err);

    _subfields = src._subfields;

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::topology::Field::~Field(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Field::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;
    err = VecDestroy(&_localVec);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_globalVec);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_outputVec);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set label for field.
void
pylith::topology::Field::setLabel(const char* value) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;

    _label = value;
    if (_dm) {
        err = PetscObjectSetName((PetscObject) _dm, value);PYLITH_CHECK_ERROR(err);
    } // of
    if (_localVec) {
        err = PetscObjectSetName((PetscObject) _localVec, value);PYLITH_CHECK_ERROR(err);
    } // if
    if (_globalVec) {
        err = PetscObjectSetName((PetscObject) _globalVec, value);PYLITH_CHECK_ERROR(err);
    } // if
    if (_outputVec) {
        err = PetscObjectSetName((PetscObject) _outputVec, value);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // setLabel


// ------------------------------------------------------------------------------------------------
// Get spatial dimension of domain.
int
pylith::topology::Field::getSpaceDim(void) const {
    const spatialdata::geocoords::CoordSys* cs = _mesh.getCoordSys();
    return (cs) ? cs->getSpaceDim() : 0;
} // getSpaceDim


// ------------------------------------------------------------------------------------------------
// Get the chart size.
PylithInt
pylith::topology::Field::chartSize(void) const {
    PYLITH_METHOD_BEGIN;

    assert(_dm);
    PetscSection s = NULL;
    PylithInt pStart, pEnd;
    PetscErrorCode err;

    err = DMGetSection(_dm, &s);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetChart(s, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(pEnd-pStart);
} // chartSize


// ------------------------------------------------------------------------------------------------
// Get the number of degrees of freedom.
PylithInt
pylith::topology::Field::getStorageSize(void) const {
    PYLITH_METHOD_BEGIN;

    PylithInt size = 0;
    if (_dm) {
        PetscSection s = NULL;
        PetscErrorCode err;
        err = DMGetSection(_dm, &s);PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetStorageSize(s, &size);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_RETURN(size);
} // getStorageSize


// ------------------------------------------------------------------------------------------------
void
pylith::topology::Field::createDiscretization(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = DMCreateDS(_dm);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // createDiscretization


// ------------------------------------------------------------------------------------------------
// Allocate PETSc section.
void
pylith::topology::Field::allocate(void) {
    PYLITH_METHOD_BEGIN;

    PetscSection s = NULL;
    PetscErrorCode err;

    /* Ideally, we should create the DS *before* adding the boundary conditions to the DS. For now, this does work.
     *
     * We cannot create the DS until after setting the discretizations for each field in subfieldsSetup() and
     * _setupLagrangeMultiplier().
     */

    err = DMGetSection(_dm, &s);PYLITH_CHECK_ERROR(err);assert(s); // Creates local section
    err = DMSetGlobalSection(_dm, NULL);PYLITH_CHECK_ERROR(err); // Creates global section
    err = PetscSectionSetUp(s);PYLITH_CHECK_ERROR(err);

    err = VecDestroy(&_localVec);PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(_dm, &_localVec);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) _localVec,  _label.c_str());PYLITH_CHECK_ERROR(err);
    err = VecSet(_localVec, 0.0);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // allocate


// ------------------------------------------------------------------------------------------------
// Zero local vector values (including constrained DOF).
void
pylith::topology::Field::zeroLocal(void) {
    PYLITH_METHOD_BEGIN;

    assert(_localVec);
    PetscErrorCode err = VecSet(_localVec, 0.0);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // zeroLocal


// ------------------------------------------------------------------------------------------------
// Print field to standard out.
void
pylith::topology::Field::view(const char* label,
                              const ViewOptions options) const {
    PYLITH_METHOD_BEGIN;

    std::cout << "Viewing field '" << _label << "' "<< label << ".\n";
    if (_subfields.size() > 0) {
        std::cout << "  Subfields:\n";
        for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
            const char* sname = s_iter->first.c_str();
            const SubfieldInfo& sinfo = s_iter->second;
            std::cout << "    Subfield " << sname << ", index: " << sinfo.index;
            const int nscomps = sinfo.description.numComponents;
            if (nscomps > 0) {
                std::cout << ", components:";
                for (int i = 0; i < nscomps; ++i) {
                    std::cout << " " << sinfo.description.componentNames[i];
                } // for
            } // if
            std::string cellBasisString = "default";
            switch (sinfo.fe.cellBasis) {
            case SIMPLEX_BASIS:
                cellBasisString = "simplex";
                break;
            case TENSOR_BASIS:
                cellBasisString = "tensor";
                break;
            case DEFAULT_BASIS:
                cellBasisString = "default";
                break;
            default:
                assert(0);
                throw std::logic_error("Unknown cell basis");
            } // switch
            std::cout << ", scale: " << sinfo.description.scale
                      << ", basisOrder: " << sinfo.fe.basisOrder
                      << ", quadOrder: " << sinfo.fe.quadOrder
                      << ", dimension: " << sinfo.fe.dimension
                      << ", cellBasis: " << cellBasisString
                      << "\n";
        } // for
    } // if

    if (_dm) {
        PetscSection section = NULL;
        PetscMPIInt numProcs, rank;
        PetscErrorCode err;

        err = DMGetSection(_dm, &section);PYLITH_CHECK_ERROR(err);
        if ((VIEW_LAYOUT == options) || (VIEW_ALL == options)) {
            err = DMView(_dm, PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
            err = PetscSectionView(section, PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
        } // if
        if ((VIEW_VALUES == options) || (VIEW_ALL == options)) {
            err = MPI_Comm_size(PetscObjectComm((PetscObject) _dm), &numProcs);PYLITH_CHECK_ERROR(err);
            err = MPI_Comm_rank(PetscObjectComm((PetscObject) _dm), &rank);PYLITH_CHECK_ERROR(err);
            for (PetscInt p = 0; p < numProcs; ++p) {
                err = PetscPrintf(PetscObjectComm((PetscObject) _dm), "Proc %d local vector\n", p);PYLITH_CHECK_ERROR(err);
                if (p == rank) {err = VecView(_localVec, PETSC_VIEWER_STDOUT_SELF);PYLITH_CHECK_ERROR(err); }
                err = PetscBarrier((PetscObject) _dm);PYLITH_CHECK_ERROR(err);
            } // for
        } // if
    } // if

    PYLITH_METHOD_END;
} // view


// ------------------------------------------------------------------------------------------------
// Scatter section information across processors to update the
// PETSc vector view of the field.
void
pylith::topology::Field::scatterLocalToVector(const PetscVec vector,
                                              InsertMode mode) const {
    PYLITH_METHOD_BEGIN;
    assert(vector);

    PetscErrorCode err;
    if (_dm) {
        assert(_localVec);
        err = DMLocalToGlobalBegin(_dm, _localVec, mode, vector);PYLITH_CHECK_ERROR(err);
        err = DMLocalToGlobalEnd(_dm, _localVec, mode, vector);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // scatterLocalToVector


// ------------------------------------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
void
pylith::topology::Field::scatterVectorToLocal(const PetscVec vector,
                                              InsertMode mode) const {
    PYLITH_METHOD_BEGIN;
    assert(vector);

    PetscErrorCode err;
    if (_dm) {
        assert(_localVec);
        err = DMGlobalToLocalBegin(_dm, vector, mode, _localVec);PYLITH_CHECK_ERROR(err);
        err = DMGlobalToLocalEnd(_dm, vector, mode, _localVec);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // scatterVectorToLocal


// ------------------------------------------------------------------------------------------------
// Scatter section information across processors to update the
// output view of the field.
void
pylith::topology::Field::scatterLocalToOutput(InsertMode mode) const {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;
    PetscDM dmOutput = NULL;
    err = DMGetOutputDM(_dm, &dmOutput);PYLITH_CHECK_ERROR(err);
    if (dmOutput) {
        assert(_localVec);
        err = DMLocalToGlobalBegin(dmOutput, _localVec, mode, _outputVec);PYLITH_CHECK_ERROR(err);
        err = DMLocalToGlobalEnd(dmOutput, _localVec, mode, _outputVec);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // scatterLocalToOutput


// ------------------------------------------------------------------------------------------------
// Add subfield to current field (inteface to use from SWIG).
void
pylith::topology::Field::subfieldAdd(const char *name,
                                     const char* alias,
                                     const VectorFieldEnum fieldType,
                                     const char* components[],
                                     const int numComponents,
                                     const double scale,
                                     const int basisOrder,
                                     const int quadOrder,
                                     const int dimension,
                                     const CellBasis cellBasis,
                                     const bool isBasisContinuous,
                                     const SpaceEnum feSpace) {
    assert(numComponents > 0);
    assert(dimension != 0);

    Description description;
    description.label = name;
    description.alias = alias;
    description.vectorFieldType = fieldType;
    description.numComponents = numComponents;
    description.componentNames.resize(numComponents);
    for (int i = 0; i < numComponents; ++i) {
        description.componentNames[i] = components[i];
    } // for
    description.scale = scale;
    description.validator = NULL;

    Discretization discretization;
    discretization.basisOrder = basisOrder;
    discretization.quadOrder = quadOrder;
    discretization.dimension = dimension;
    discretization.cellBasis = cellBasis;
    discretization.isBasisContinuous = isBasisContinuous;
    discretization.feSpace = feSpace;

    this->subfieldAdd(description, discretization);
} // subfieldAdd


// ------------------------------------------------------------------------------------------------
// Add subfield.
void
pylith::topology::Field::subfieldAdd(const Description& description,
                                     const Discretization& discretization) {
    PYLITH_METHOD_BEGIN;

    assert(0 == _subfields.count(description.label));

    // Keep track of name/components until setup
    SubfieldInfo info;
    info.description = description;
    info.fe = discretization;
    info.index = _subfields.size(); // Indices match order added.
    _subfields[description.label] = info;

    PYLITH_METHOD_END;
} // subfieldAdd


// ------------------------------------------------------------------------------------------------
void
pylith::topology::Field::subfieldsSetup(void) {
    PYLITH_METHOD_BEGIN;

    assert(_dm);

    // Setup section now that we know the total number of sub-fields and components.
    PetscErrorCode err;

    bool quadOrderSet = false;
    int quadOrder = -999;
    for (subfields_type::iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
        const char* sname = s_iter->first.c_str();
        SubfieldInfo& sinfo = s_iter->second;

        if (quadOrderSet) {
            if (quadOrder != sinfo.fe.quadOrder) {
                std::ostringstream msg;
                msg << "PETSc DMPlex routines currently assume all subfields use the same quadrature order. Quaadrature order of "
                    << sinfo.fe.quadOrder << " for subfield '" << sname << "' does not match the quadrature order of " << quadOrder
                    << " for other subfields in field '" << getLabel() << "'.";
                throw std::runtime_error(msg.str());
            } // if
        } else {
            quadOrder = sinfo.fe.quadOrder;
            quadOrderSet = true;
        } // if/else

        sinfo.fe.cellBasis = _mesh.isSimplex() ? SIMPLEX_BASIS : TENSOR_BASIS;
        PetscFE fe = FieldOps::createFE(sinfo.fe, _dm, sinfo.description.numComponents);assert(fe);
        err = PetscFESetName(fe, sname);PYLITH_CHECK_ERROR(err);
        err = DMSetField(_dm, sinfo.index, NULL, (PetscObject)fe);PYLITH_CHECK_ERROR(err);
        err = DMSetFieldAvoidTensor(_dm, sinfo.index, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
        err = PetscFEDestroy(&fe);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // subfieldsSetup


// ------------------------------------------------------------------------------------------------
// Does field have given auxiliary subfield?
bool
pylith::topology::Field::hasSubfield(const char* name) const {
    PYLITH_METHOD_BEGIN;

    subfields_type::const_iterator iter = _subfields.find(name);
    PYLITH_METHOD_RETURN(_subfields.end() != iter);
} // hasSubfield


// ------------------------------------------------------------------------------------------------
// Get names of subfields.
pylith::string_vector
pylith::topology::Field::subfieldNames(void) const {
    PYLITH_METHOD_BEGIN;

    const size_t numSubfields = _subfields.size();
    pylith::string_vector names(numSubfields);
    for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
        const SubfieldInfo& sinfo = s_iter->second;
        names[sinfo.index] = s_iter->first;
    } // for

    PYLITH_METHOD_RETURN(pylith::string_vector(names));
} // subfieldNames


// ------------------------------------------------------------------------------------------------
// Get metadata for subfield.
const pylith::topology::Field::SubfieldInfo&
pylith::topology::Field::subfieldInfo(const char* name) const {
    PYLITH_METHOD_BEGIN;

    subfields_type::const_iterator iter = _subfields.find(name);
    if (_subfields.end() == iter) {
        std::ostringstream msg;
        msg << "Could not find subfield '" << name << "' in field '" << getLabel() << "'.\n"
            << "Available subfields:";
        for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
            msg << " '" << s_iter->first << "'";
        } // for
        msg << std::endl;

        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(iter->second);
} // subfieldInfo


// ------------------------------------------------------------------------------------------------
// Create global vector.
void
pylith::topology::Field::createGlobalVector(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = VecDestroy(&_globalVec);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(_dm, &_globalVec);PYLITH_CHECK_ERROR(err);assert(_globalVec);
    err = PetscObjectSetName((PetscObject) _globalVec, getLabel());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Create global vector with no constrained DOF.
void
pylith::topology::Field::createOutputVector(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = VecDestroy(&_outputVec);PYLITH_CHECK_ERROR(err);

    PetscDM dmOutput = NULL;
    err = DMGetOutputDM(_dm, &dmOutput);PYLITH_CHECK_ERROR(err);
    PetscDS dsOutput = NULL;
    err = DMGetDS(dmOutput, &dsOutput);PYLITH_CHECK_ERROR(err);
    err = PetscDSSetUp(dsOutput);PYLITH_CHECK_ERROR(err);

    err = DMCreateGlobalVector(dmOutput, &_outputVec);PYLITH_CHECK_ERROR(err);assert(_outputVec);
    err = PetscObjectSetName((PetscObject) _outputVec, getLabel());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
}


// End of file
