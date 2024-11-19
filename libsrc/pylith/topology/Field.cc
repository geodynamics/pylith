// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/topology/Field.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/topology/MeshOps.hh" // USES isCohesiveCell()
#include "pylith/faults/TopologyOps.hh" // USES getInterfacesLabel()

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include <cassert> // USES assert()
#include <iostream> // USES std::cout
#include <iomanip> // USES std::setw()

namespace pylith {
    namespace topology {
        class _Field {
public:

            /// View field layout.
            static
            void viewLayout(const pylith::topology::Field& field);

            /// View values in local vector of field.
            static
            void viewValues(const pylith::topology::Field& field);

            /// View local section.
            static
            void viewLocalSection(const Field& field);

            /// View global section.
            static
            void viewGlobalSection(const Field& field);

        }; // _Field
    }
}

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::topology::Field::Field(const pylith::topology::Mesh& mesh) :
    _mesh(NULL),
    _localVec(NULL),
    _globalVec(NULL),
    _outputVec(NULL) {
    PYLITH_METHOD_BEGIN;

    GenericComponent::setName("field");
    _label = "unknown";
    _mesh = mesh.clone();assert(_mesh);

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Constructor with field to use for layout.
pylith::topology::Field::Field(const Field& src) :
    _mesh(NULL),
    _localVec(NULL),
    _globalVec(NULL),
    _outputVec(NULL) {
    PYLITH_METHOD_BEGIN;

    _subfields = src._subfields;

    if (!src._mesh) {
        PYLITH_JOURNAL_LOGICERROR("Source field _mesh must be non-NULL.");
    } // if
    _mesh = src._mesh->clone();

    PetscErrorCode err;
    assert(_mesh);
    err = DMCopyDisc(src._mesh->getDM(), _mesh->getDM());PYLITH_CHECK_ERROR(err);

    assert(!_localVec);
    err = DMCreateLocalVector(_mesh->getDM(), &_localVec);PYLITH_CHECK_ERROR(err);
    err = VecSet(_localVec, 0.0);PYLITH_CHECK_ERROR(err);

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

    delete _mesh;_mesh = NULL;

    PetscErrorCode err;
    err = VecDestroy(&_localVec);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_globalVec);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_outputVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Get mesh associated with field.
const pylith::topology::Mesh&
pylith::topology::Field::getMesh(void) const {
    assert(_mesh);
    return *_mesh;
}


// ------------------------------------------------------------------------------------------------
// Get PETSc DM associated with field.
PetscDM
pylith::topology::Field::getDM(void) const {
    return (_mesh) ? _mesh->getDM() : NULL;
}


// ------------------------------------------------------------------------------------------------
// Get label for field.
const char*
pylith::topology::Field::getLabel(void) const {
    return _label.c_str();
}


// ------------------------------------------------------------------------------------------------
// Set label for field.
void
pylith::topology::Field::setLabel(const char* value) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PetscErrorCode err;

    _label = value;
    if (_mesh->getDM()) {
        err = PetscObjectSetName((PetscObject) _mesh->getDM(), value);PYLITH_CHECK_ERROR(err);
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
// Get local PetscSection.
PetscSection
pylith::topology::Field::getLocalSection(void) const {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PetscSection s = NULL;
    PetscErrorCode err = DMGetSection(_mesh->getDM(), &s);PYLITH_CHECK_ERROR(err);
    PYLITH_METHOD_RETURN(s);
}


// ------------------------------------------------------------------------------------------------
// Get global PetscSection.
PetscSection
pylith::topology::Field::getGlobalSection(void) const {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PetscSection s = NULL;
    PetscErrorCode err = DMGetGlobalSection(_mesh->getDM(), &s);PYLITH_CHECK_ERROR(err);
    PYLITH_METHOD_RETURN(s);
}


// ------------------------------------------------------------------------------------------------
// Get local vector.
PetscVec
pylith::topology::Field::getLocalVector(void) const {
    return _localVec;
}


// ------------------------------------------------------------------------------------------------
// Get global vector.
PetscVec
pylith::topology::Field::getGlobalVector(void) const {
    return _globalVec;
}


// ------------------------------------------------------------------------------------------------
// Get global vector without constrained degrees of freedom for output.
PetscVec
pylith::topology::Field::getOutputVector(void) const {
    return _outputVec;
}


// ------------------------------------------------------------------------------------------------
// Get spatial dimension of coordinate system for field.
size_t
pylith::topology::Field::getSpaceDim(void) const {
    assert(_mesh);
    const spatialdata::geocoords::CoordSys* cs = _mesh->getCoordSys();assert(cs);
    return cs->getSpaceDim();
}


// ------------------------------------------------------------------------------------------------
// Get the chart size.
PylithInt
pylith::topology::Field::getChartSize(void) const {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PetscSection s = NULL;
    PylithInt pStart, pEnd;
    PetscErrorCode err;

    err = DMGetSection(_mesh->getDM(), &s);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetChart(s, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(pEnd-pStart);
} // chartSize


// ------------------------------------------------------------------------------------------------
// Get the number of degrees of freedom.
PylithInt
pylith::topology::Field::getStorageSize(void) const {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PylithInt size = 0;
    PetscSection s = NULL;
    PetscErrorCode err;
    err = DMGetSection(_mesh->getDM(), &s);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetStorageSize(s, &size);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(size);
} // getStorageSize


// ------------------------------------------------------------------------------------------------
void
pylith::topology::Field::createDiscretization(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PetscErrorCode err = DMCreateDS(_mesh->getDM());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // createDiscretization


// ------------------------------------------------------------------------------------------------
// Allocate PETSc section.
void
pylith::topology::Field::allocate(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PetscSection s = NULL;
    PetscErrorCode err;

    /* Ideally, we should create the DS *before* adding the boundary conditions to the DS. For now, this does work.
     *
     * We cannot create the DS until after setting the discretizations for each field in subfieldsSetup() and
     * _setupLagrangeMultiplier().
     */

    const PetscDM dm = _mesh->getDM();assert(dm);
    err = DMGetSection(dm, &s);PYLITH_CHECK_ERROR(err);assert(s); // Creates local section
    err = DMSetGlobalSection(dm, NULL);PYLITH_CHECK_ERROR(err); // Creates global section
    err = PetscSectionSetUp(s);PYLITH_CHECK_ERROR(err);

    err = VecDestroy(&_localVec);PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(dm, &_localVec);PYLITH_CHECK_ERROR(err);
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
    assert(_mesh);

    const int rank = _mesh->getCommRank();
    if (!rank) {
        std::cout << "Viewing field '" << _label << "' "<< label << ".\n";
        if (_subfields.size() > 0) {
            std::cout << "  Subfields:\n";
            for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
                const char* sname = s_iter->first.c_str();
                const SubfieldInfo& sinfo = s_iter->second;
                std::cout << "    Subfield " << sname << ", index: " << sinfo.index;
                const int numComponents = sinfo.description.numComponents;
                if (numComponents > 0) {
                    std::cout << ", components:";
                    for (int i = 0; i < numComponents; ++i) {
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
                    PYLITH_JOURNAL_LOGICERROR("Unknown cell basis");
                } // switch
                std::cout << ", scale: " << sinfo.description.scale
                          << ", basisOrder: " << sinfo.fe.basisOrder
                          << ", quadOrder: " << sinfo.fe.quadOrder
                          << ", dimension: " << sinfo.fe.dimension
                          << ", cellBasis: " << cellBasisString
                          << "\n";
            } // for
        } // if
    } // if

    if ((VIEW_LAYOUT == options) || (VIEW_ALL == options)) {
        _Field::viewLayout(*this);
    } // if
    const PetscDM dm = _mesh->getDM();assert(dm);
    PetscErrorCode err = PetscBarrier((PetscObject) dm);PYLITH_CHECK_ERROR(err);

    if ((VIEW_VALUES == options) || (VIEW_ALL == options)) {
        _Field::viewValues(*this);
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
    assert(_mesh);
    assert(vector);

    PetscErrorCode err;
    assert(_localVec);
    err = DMLocalToGlobalBegin(_mesh->getDM(), _localVec, mode, vector);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(_mesh->getDM(), _localVec, mode, vector);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // scatterLocalToVector


// ------------------------------------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
void
pylith::topology::Field::scatterVectorToLocal(const PetscVec vector,
                                              InsertMode mode) const {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(vector);

    PetscErrorCode err;
    assert(_localVec);
    err = DMGlobalToLocalBegin(_mesh->getDM(), vector, mode, _localVec);PYLITH_CHECK_ERROR(err);
    err = DMGlobalToLocalEnd(_mesh->getDM(), vector, mode, _localVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // scatterVectorToLocal


// ------------------------------------------------------------------------------------------------
// Scatter section information across processors to update the
// output view of the field.
void
pylith::topology::Field::scatterLocalToOutput(InsertMode mode) const {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PetscErrorCode err;
    PetscDM dmOutput = NULL;
    err = DMGetOutputDM(_mesh->getDM(), &dmOutput);PYLITH_CHECK_ERROR(err);
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
                                     const bool isFaultOnly,
                                     const CellBasis cellBasis,
                                     const SpaceEnum feSpace,
                                     const bool isBasisContinuous) {
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
    discretization.feSpace = feSpace;
    discretization.isBasisContinuous = isBasisContinuous;
    discretization.isFaultOnly = isFaultOnly;

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
    assert(_mesh);

    // Setup section now that we know the total number of sub-fields and components.
    PetscErrorCode err;

    CellBasis cellBasis = pylith::topology::MeshOps::isSimplexMesh(*_mesh) ? SIMPLEX_BASIS : TENSOR_BASIS;
    const PetscDM dm = _mesh->getDM();

    bool quadOrderSet = false;
    int quadOrder = -999;
    for (subfields_type::iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
        const char* sname = s_iter->first.c_str();
        SubfieldInfo& sinfo = s_iter->second;

        if (quadOrderSet) {
            if (quadOrder != sinfo.fe.quadOrder) {
                std::ostringstream msg;
                msg << "PETSc DMPlex routines currently assume all subfields use the same quadrature order. Quadrature order of "
                    << sinfo.fe.quadOrder << " for subfield '" << sname << "' does not match the quadrature order of " << quadOrder
                    << " for other subfields in field '" << getLabel() << "'.";
                throw std::runtime_error(msg.str());
            } // if
        } else {
            quadOrder = sinfo.fe.quadOrder;
            quadOrderSet = true;
        } // if/else

        sinfo.fe.cellBasis = cellBasis;
        PetscFE fe = FieldOps::createFE(sinfo.fe, dm, sinfo.description.numComponents);assert(fe);
        err = PetscFESetName(fe, sname);PYLITH_CHECK_ERROR(err);

        // :KLUDGE: We need PETSc DMPlex to support specifying subfields over labels and label values.
        // Once that is implemented, then we should switch to specifying subfields over labels and values.
        // For now we assume subfields are on:
        //   + all degrees of freedom,
        //   + everywhere but fault degrees of freedom, or
        //   + only fault degrees of freedom.
        if (!sinfo.fe.isFaultOnly) {
            err = DMSetField(dm, sinfo.index, NULL, (PetscObject)fe);PYLITH_CHECK_ERROR(err);
            err = DMSetFieldAvoidTensor(dm, sinfo.index, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
        } else {
            PetscDMLabel interfacesLabel = pylith::faults::TopologyOps::getInterfacesLabel(dm);
            err = DMSetField(dm, sinfo.index, interfacesLabel, (PetscObject)fe);PYLITH_CHECK_ERROR(err);
        } // if/else
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
pylith::topology::Field::getSubfieldNames(void) const {
    PYLITH_METHOD_BEGIN;

    const size_t numSubfields = _subfields.size();
    pylith::string_vector names(numSubfields);
    for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
        const SubfieldInfo& sinfo = s_iter->second;
        names[sinfo.index] = s_iter->first;
    } // for

    PYLITH_METHOD_RETURN(names);
} // subfieldNames


// ------------------------------------------------------------------------------------------------
// Get metadata for subfield.
const pylith::topology::Field::SubfieldInfo&
pylith::topology::Field::getSubfieldInfo(const char* name) const {
    PYLITH_METHOD_BEGIN;
    assert(name);

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
    assert(_mesh);

    PetscErrorCode err = VecDestroy(&_globalVec);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(_mesh->getDM(), &_globalVec);PYLITH_CHECK_ERROR(err);assert(_globalVec);
    err = PetscObjectSetName((PetscObject) _globalVec, getLabel());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Create global vector with no constrained DOF.
void
pylith::topology::Field::createOutputVector(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PetscErrorCode err = VecDestroy(&_outputVec);PYLITH_CHECK_ERROR(err);

    PetscDM dmOutput = NULL;
    err = DMGetOutputDM(_mesh->getDM(), &dmOutput);PYLITH_CHECK_ERROR(err);

    PetscDS dsOutput = NULL;
    PetscInt numRegions = 0;
    err = DMGetNumDS(dmOutput, &numRegions);PYLITH_CHECK_ERROR(err);
    for (PetscInt i = 0; i < numRegions; ++i) {
        err = DMGetRegionNumDS(dmOutput, i, NULL, NULL, &dsOutput, NULL);PYLITH_CHECK_ERROR(err);
        err = PetscDSSetUp(dsOutput);PYLITH_CHECK_ERROR(err);
    } // for

    err = DMCreateGlobalVector(dmOutput, &_outputVec);PYLITH_CHECK_ERROR(err);assert(_outputVec);
    err = PetscObjectSetName((PetscObject) _outputVec, getLabel());PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// View field layout.
void
pylith::topology::_Field::viewLayout(const Field& field) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = PETSC_SUCCESS;
    const PetscDM dm = field.getDM();assert(dm);

    PetscMPIInt numProcs, rank;
    err = MPI_Comm_size(PetscObjectComm((PetscObject) dm), &numProcs);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(PetscObjectComm((PetscObject) dm), &rank);PYLITH_CHECK_ERROR(err);

    const pylith::string_vector& subfieldNames = field.getSubfieldNames();
    const size_t numSubfields = subfieldNames.size();
    if (!rank) {
        std::cout << "      ";
        for (size_t i = 0; i < numSubfields; ++i) {
            std::cout << "  " << subfieldNames[i];
        } // for`
        std::cout << "\n";
    } // if

    _Field::viewLocalSection(field);
    _Field::viewGlobalSection(field);

    PYLITH_METHOD_END;
} // viewLayout


// ------------------------------------------------------------------------------------------------
// View values in local vector of field.
void
pylith::topology::_Field::viewValues(const Field& field) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = PETSC_SUCCESS;
    const PetscDM dm = field.getDM();assert(dm);

    PetscMPIInt numProcs, rank;
    err = MPI_Comm_size(PetscObjectComm((PetscObject) dm), &numProcs);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(PetscObjectComm((PetscObject) dm), &rank);PYLITH_CHECK_ERROR(err);

    const pylith::string_vector& subfieldNames = field.getSubfieldNames();
    const size_t numSubfields = subfieldNames.size();
    pylith::int_array numComponents(numSubfields);
    for (size_t i = 0; i < numSubfields; ++i) {
        numComponents[i] = field.getSubfieldInfo(subfieldNames[i].c_str()).description.numComponents;
    } // for

    const size_t pwidth = 6;
    const size_t fwidth = 12;
    std::cout << std::scientific << std::setprecision(4);
    const char* none = "--";
    VecVisitorMesh fieldVisitor(field);
    PetscInt pStart, pEnd;
    err = PetscSectionGetChart(fieldVisitor.selectedSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    const PetscScalar* localArray = fieldVisitor.localArray();
    for (PetscInt iProc = 0; iProc < numProcs; ++iProc) {
        std::cout << "Processor " << iProc << "\n";
        for (PetscInt point = pStart; point < pEnd; ++point) {
            std::cout << std::setw(pwidth) << point;
            for (size_t iField = 0; iField < numSubfields; ++iField) {
                const PetscInt numDof = fieldVisitor.sectionSubfieldDof(iField, point);
                if (numDof > 0) {
                    assert(numDof == numComponents[iField]);
                    const PetscInt off = fieldVisitor.sectionSubfieldOffset(iField, point);
                    for (PetscInt iDof = 0; iDof < numDof; ++iDof) {
                        std::cout << std::setw(fwidth) << localArray[off+iDof];
                    } // for
                } else {
                    for (PetscInt iDof = 0; iDof < numComponents[iField]; ++iDof) {
                        std::cout << std::setw(fwidth) << none;
                    } // for
                } // if/else
                std::cout << "  ";
            } // for
            std::cout << "\n";
        } // for
        err = PetscBarrier((PetscObject) dm);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// View local section layout.
void
pylith::topology::_Field::viewLocalSection(const Field& field) {
    PYLITH_METHOD_BEGIN;
    const size_t pwidth = 6;
    const char* none = "--";

    PetscErrorCode err = PETSC_SUCCESS;
    const PetscDM dm = field.getDM();assert(dm);
    PetscMPIInt numProcs, rank;
    err = MPI_Comm_size(PetscObjectComm((PetscObject) dm), &numProcs);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(PetscObjectComm((PetscObject) dm), &rank);PYLITH_CHECK_ERROR(err);

    const pylith::string_vector& subfieldNames = field.getSubfieldNames();
    const size_t numSubfields = subfieldNames.size();
    pylith::int_array numComponents(numSubfields);
    for (size_t i = 0; i < numSubfields; ++i) {
        numComponents[i] = field.getSubfieldInfo(subfieldNames[i].c_str()).description.numComponents;
    } // for

    PetscInt pStart, pEnd;
    pylith::topology::VecVisitorMesh fieldVisitor(field);
    fieldVisitor.selectSection(pylith::topology::VecVisitorMesh::LOCAL_SECTION);
    const PetscSection section = fieldVisitor.selectedSection();
    err = PetscSectionGetChart(section, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    if (!rank) {
        std::cout << "Local Section" << "\n";
    } // if

    for (PetscInt iProc = 0; iProc < numProcs; ++iProc) {
        std::cout << "Processor " << iProc << "\n";
        for (PetscInt point = pStart; point < pEnd; ++point) {
            std::cout << std::setw(pwidth) << point;
            for (size_t iField = 0; iField < numSubfields; ++iField) {
                const PetscInt numDof = fieldVisitor.sectionSubfieldDof(iField, point);
                if (numDof > 0) {
                    const PetscInt off = fieldVisitor.sectionSubfieldOffset(iField, point);
                    for (PetscInt iDof = 0; iDof < numDof; ++iDof) {
                        std::cout << std::setw(pwidth) << off+iDof;
                    } // for
                } else {
                    for (PetscInt iDof = 0; iDof < numComponents[iField]; ++iDof) {
                        std::cout << std::setw(pwidth) << none;
                    } // for
                } // if/else
                std::cout << "  ";
            } // for
            PetscInt numConstraintDof;
            err = PetscSectionGetConstraintDof(section, point, &numConstraintDof);PYLITH_CHECK_ERROR(err);
            if (numConstraintDof > 0) {
                const PetscInt off = fieldVisitor.sectionOffset(point);
                std::cout << "  (constrained:";
                const PetscInt* constraintIndices = NULL;
                err = PetscSectionGetConstraintIndices(section, point, &constraintIndices);PYLITH_CHECK_ERROR(err);
                for (PetscInt iDof = 0; iDof < numConstraintDof; ++iDof) {
                    std::cout << " " << off+constraintIndices[iDof];
                } // for
                std::cout << ")";
            }
            std::cout << "\n";
        } // for
        err = PetscBarrier((PetscObject) dm);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // _viewLocalLayout


// ------------------------------------------------------------------------------------------------
// View global section layout.
void
pylith::topology::_Field::viewGlobalSection(const Field& field) {
    PYLITH_METHOD_BEGIN;
    const size_t pwidth = 6;
    const char* none = "--";

    PetscErrorCode err = PETSC_SUCCESS;
    const PetscDM dm = field.getDM();assert(dm);
    PetscMPIInt numProcs, rank;
    err = MPI_Comm_size(PetscObjectComm((PetscObject) dm), &numProcs);PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(PetscObjectComm((PetscObject) dm), &rank);PYLITH_CHECK_ERROR(err);

    const pylith::string_vector& subfieldNames = field.getSubfieldNames();
    const size_t numSubfields = subfieldNames.size();
    pylith::int_array numComponents(numSubfields);
    for (size_t i = 0; i < numSubfields; ++i) {
        numComponents[i] = field.getSubfieldInfo(subfieldNames[i].c_str()).description.numComponents;
    } // for

    PetscInt pStart, pEnd;
    pylith::topology::VecVisitorMesh fieldVisitor(field);
    fieldVisitor.selectSection(pylith::topology::VecVisitorMesh::GLOBAL_SECTION);
    const PetscSection section = fieldVisitor.selectedSection();
    err = PetscSectionGetChart(section, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    if (!rank) {
        std::cout << "Global Section" << "\n";
    } // if

    for (PetscInt iProc = 0; iProc < numProcs; ++iProc) {
        std::cout << "Processor " << iProc << "\n";
        for (PetscInt point = pStart; point < pEnd; ++point) {
            std::cout << std::setw(pwidth) << point;
            for (size_t iField = 0; iField < numSubfields; ++iField) {
                const PetscInt numDof = fieldVisitor.sectionSubfieldDof(iField, point);
                if (numDof > 0) {
                    PetscInt numConstraintDof;
                    err = PetscSectionGetFieldConstraintDof(section, point, iField, &numConstraintDof);PYLITH_CHECK_ERROR(err);
                    if (numConstraintDof > 0) {
                        assert(numConstraintDof <= numDof);
                        PetscInt iConstraint = 0;
                        const PetscInt* constraintIndices = NULL;
                        err = PetscSectionGetFieldConstraintIndices(section, point, iField, &constraintIndices);PYLITH_CHECK_ERROR(err);
                        const PetscInt off = fieldVisitor.sectionSubfieldOffset(iField, point);
                        for (PetscInt iDof = 0; iDof < numDof; ++iDof) {
                            if ((iConstraint < numConstraintDof) && (constraintIndices[iConstraint] == iDof)) {
                                std::cout << std::setw(pwidth) << none;
                                ++iConstraint;
                            } else {
                                std::cout << std::setw(pwidth) << off+iDof-iConstraint;
                            } // if/else
                        } // for
                    } else {
                        const PetscInt off = fieldVisitor.sectionSubfieldOffset(iField, point);
                        for (PetscInt iDof = 0; iDof < numDof; ++iDof) {
                            std::cout << std::setw(pwidth) << off+iDof;
                        } // for
                    } // if/else
                } else {
                    for (PetscInt iDof = 0; iDof < numComponents[iField]; ++iDof) {
                        std::cout << std::setw(pwidth) << none;
                    } // for
                } // if/else
                std::cout << "  ";
            } // for
            std::cout << "\n";
        } // for
        err = PetscBarrier((PetscObject) dm);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // _viewGlobalLayout


// End of file
