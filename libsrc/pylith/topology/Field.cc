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

#include "petscds.h" // USES PetscDS

#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include <iostream>

// ----------------------------------------------------------------------
const char* pylith::topology::Field::_genericComponent = "field";

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Field::Field(const Mesh& mesh) :
    _mesh(mesh),
    _dm(NULL),
    _localVec(NULL)
{ // constructor
    PYLITH_METHOD_BEGIN;

    GenericComponent::name(_genericComponent);

    _dimsOkay = false;

    _label = "unknown";
    if (mesh.dmMesh()) {
        PetscDM dm = mesh.dmMesh(); assert(dm);
        PetscVec coordVec = NULL;
        PetscSection s = NULL;
        PetscErrorCode err;

        err = DMDestroy(&_dm); PYLITH_CHECK_ERROR(err);
        err = DMClone(dm, &_dm); PYLITH_CHECK_ERROR(err);
        err = DMGetCoordinatesLocal(dm, &coordVec); PYLITH_CHECK_ERROR(err);
        if (coordVec) {
            PetscDM coordDM = NULL, newCoordDM = NULL;
            PetscSection coordSection = NULL, newCoordSection = NULL;

            err = DMGetCoordinateDM(dm, &coordDM); PYLITH_CHECK_ERROR(err);
            err = DMGetCoordinateDM(_dm, &newCoordDM); PYLITH_CHECK_ERROR(err);
            err = DMGetDefaultSection(coordDM, &coordSection); PYLITH_CHECK_ERROR(err);
            err = PetscSectionClone(coordSection, &newCoordSection); PYLITH_CHECK_ERROR(err);
            err = DMSetDefaultSection(newCoordDM, newCoordSection); PYLITH_CHECK_ERROR(err);
            err = PetscSectionDestroy(&newCoordSection); PYLITH_CHECK_ERROR(err);
            err = DMSetCoordinatesLocal(_dm, coordVec); PYLITH_CHECK_ERROR(err);
        } // if

#if 1 // :TEMPORARY: These 3 lines go away once we have everything converted to using PetscDS.
        err = PetscSectionCreate(mesh.comm(), &s); PYLITH_CHECK_ERROR(err);
        err = DMSetDefaultSection(_dm, s); PYLITH_CHECK_ERROR(err);
        err = PetscSectionDestroy(&s); PYLITH_CHECK_ERROR(err);
#endif
    } // if

    PYLITH_METHOD_END;
} // constructor

#if 0
// ----------------------------------------------------------------------
// Constructor with mesh, DM, and metadata
pylith::topology::Field::Field(const Mesh& mesh,
                               PetscDM dm,
                               const char* label) :
    _mesh(mesh),
    _dm(dm),
    _localVec(NULL)
{ // constructor
    PYLITH_METHOD_BEGIN;

    GenericComponent::name(_genericComponent);

    assert(dm);
    PetscErrorCode err;
    err = DMCreateLocalVector(_dm, &_localVec); PYLITH_CHECK_ERROR(err);

    assert(label);
    this->label(label);

    PYLITH_METHOD_END;
} // constructor
#endif

// ----------------------------------------------------------------------
// Constructor with mesh, DM, local data, and metadata
pylith::topology::Field::Field(const Mesh& mesh,
                               PetscDM dm,
                               PetscVec localVec,
                               const Description& description) :
    _label(description.label),
    _mesh(mesh),
    _dm(dm),
    _localVec(NULL)
{ // constructor
    PYLITH_METHOD_BEGIN;

    GenericComponent::name(_genericComponent);

    assert(dm);
    assert(localVec);

    PetscErrorCode err;

    Discretization fe;
    fe.basisOrder = 0;
    fe.quadOrder = 0;
    subfieldAdd(description, fe);
    subfieldsSetup();

    err = DMCreateLocalVector(_dm, &_localVec); PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) _localVec,  _label.c_str()); PYLITH_CHECK_ERROR(err);
    err = VecCopy(localVec, _localVec); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Field::~Field(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Field::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    clear();
    PetscErrorCode err = DMDestroy(&_dm); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set label for field.
void
pylith::topology::Field::label(const char* value)
{ // label
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;

    _label = value;
    if (_localVec) {
        err = PetscObjectSetName((PetscObject) _localVec, value); PYLITH_CHECK_ERROR(err);
    } // if

    const scatter_map_type::const_iterator scattersEnd = _scatters.end();
    for (scatter_map_type::const_iterator s_iter = _scatters.begin();
         s_iter != scattersEnd;
         ++s_iter) {
        if (s_iter->second.vector) {
            err = PetscObjectSetName((PetscObject)s_iter->second.vector, value); PYLITH_CHECK_ERROR(err);
        } // if
    } // for

    PYLITH_METHOD_END;
} // label

// ----------------------------------------------------------------------
// Get spatial dimension of domain.
int
pylith::topology::Field::spaceDim(void) const
{ // spaceDim
    const spatialdata::geocoords::CoordSys* cs = _mesh.coordsys();
    return (cs) ? cs->spaceDim() : 0;
} // spaceDim

// ----------------------------------------------------------------------
// Get vector field type.
pylith::topology::Field::VectorFieldEnum
pylith::topology::Field::vectorFieldType(void) const
{ // vectorFieldType
    VectorFieldEnum value = OTHER;
    if (1 == _subfields.size()) {
        value = _subfields.begin()->second.description.vectorFieldType;
    } // if

    return value;
} // vectorFieldType

// ----------------------------------------------------------------------
// Has section been setup?
bool
pylith::topology::Field::hasSection(void) const
{ // hasSection
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;
    PetscSection s = NULL;
    err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);

    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(s, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);

    bool result = (pEnd < 0) ? false : true;

    PYLITH_METHOD_RETURN(result);
} // hasSection

// ----------------------------------------------------------------------
// Get the chart size.
int
pylith::topology::Field::chartSize(void) const
{ // chartSize
    PYLITH_METHOD_BEGIN;

    assert(_dm);
    PetscSection s = NULL;
    PetscInt pStart, pEnd;
    PetscErrorCode err;

    err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetChart(s, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(pEnd-pStart);
} // chartSize

// ----------------------------------------------------------------------
// Get the number of degrees of freedom.
int
pylith::topology::Field::sectionSize(void) const
{ // sectionSize
    PYLITH_METHOD_BEGIN;

    PetscInt size = 0;

    if (_dm) {
        PetscSection s = NULL;
        PetscErrorCode err;

        err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetStorageSize(s, &size); PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_RETURN(size);
} // sectionSize

// ----------------------------------------------------------------------
// Set chart for solution.
void
pylith::topology::Field::setupSolnChart(void)
{ // setupSolnChart
    PYLITH_METHOD_BEGIN;

    assert(_dm);

    // :KLUDGE: Assume solution has DOF over vertices and hybrid edges.
    PetscErrorCode err;
    // Get range of vertices.
    PetscInt pStart = -1;
    PetscInt pEnd = -1;
    err = DMPlexGetDepthStratum(_dm, 0, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
    // Get last edge and
    PetscInt eEnd = -1;
    PetscInt eMax = -1;
    err = DMPlexGetDepthStratum(_dm, 1, NULL, &eEnd); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHybridBounds(_dm, NULL, NULL, &eMax, NULL); PYLITH_CHECK_ERROR(err);
    // If have hybrid edges, extend chart to include hybrid edges; otherwise just use points.
    if (eEnd > eMax) {
        pEnd = eEnd;
    } // if

    PetscSection s = NULL;
    if (pStart < pEnd) {
        err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
        err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetChart(s, pStart, pEnd); PYLITH_CHECK_ERROR(err);
    } else { // create empty chart
        err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
        err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetChart(s, 0, 0); PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // setupSolnChart


// ----------------------------------------------------------------------
// Setup default DOF for solution.
void
pylith::topology::Field::setupSolnDof(const int fiberDim,
                                      const char* subfieldName)
{ // setupSolnDof
    PYLITH_METHOD_BEGIN;

    assert(_dm);

    // :TODO: Update this to use discretization information after removing FIAT.

    // :KLUDGE: Assume solution has DOF over vertices.
    const int subfieldIndex = _subfields[subfieldName].index;

    PetscErrorCode err;
    // Get range of vertices.
    PetscInt pStart = -1;
    PetscInt pEnd = -1;
    err = DMPlexGetDepthStratum(_dm, 0, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);

    PetscSection s = NULL;
    err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
    for (PetscInt p = pStart; p < pEnd; ++p) {
        err = PetscSectionSetDof(s, p, fiberDim); PYLITH_CHECK_ERROR(err);

        // Set DOF in subfield
        err = PetscSectionSetFieldDof(s, p, subfieldIndex, fiberDim); PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // setupSolnDof


// ----------------------------------------------------------------------
// Create PETSc section and set chart and fiber dimesion for a list of
// points.
void
pylith::topology::Field::newSection(const int_array& points,
                                    const int fiberDim)
{ // newSection
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;

    // Clear memory
    clear();
    assert(_dm);
    if (fiberDim < 0) {
        std::ostringstream msg;
        msg << "Fiber dimension (" << fiberDim << ") for field '" << _label
            << "' must be nonnegative.";
        throw std::runtime_error(msg.str());
    } // if

    const PetscInt npts = points.size();
    if (npts > 0) {
        PetscSection s = NULL;
        PetscInt pointMin = points[0], pointMax = points[0];

        for (PetscInt i = 1; i < npts; ++i) {
            pointMin = std::min(pointMin, points[i]);
            pointMax = std::max(pointMax, points[i]);
        } // for
        err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
        err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetChart(s, pointMin, pointMax+1); PYLITH_CHECK_ERROR(err);
        for (PetscInt i = 0; i < npts; ++i) {
            err = PetscSectionSetDof(s, points[i], fiberDim); PYLITH_CHECK_ERROR(err);
        } // for
    } else { // create empty chart
        PetscSection s = NULL;

        err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
        err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetChart(s, 0, 0); PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create PETSc section and set chart and fiber dimesion for a list of
// points.
void
pylith::topology::Field::newSection(const PetscInt *points,
                                    const PetscInt num,
                                    const int fiberDim)
{ // newSection
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;

    // Clear memory
    clear();
    assert(_dm);
    if (fiberDim < 0) {
        std::ostringstream msg;
        msg << "Fiber dimension (" << fiberDim << ") for field '" << _label
            << "' must be nonnegative.";
        throw std::runtime_error(msg.str());
    } // if

    if (num > 0) {
        PetscSection s = NULL;
        PetscInt pointMin = points[0], pointMax = points[0];

        for (PetscInt i = 1; i < num; ++i) {
            pointMin = std::min(pointMin, points[i]);
            pointMax = std::max(pointMax, points[i]);
        } // for
        err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
        err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetChart(s, pointMin, pointMax+1); PYLITH_CHECK_ERROR(err);
        for (PetscInt i = 0; i < num; ++i) {
            err = PetscSectionSetDof(s, points[i], fiberDim); PYLITH_CHECK_ERROR(err);
        } // for
    } else { // create empty chart
        PetscSection s = NULL;

        err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
        err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetChart(s, 0, 0); PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create PETSc section and set chart and fiber dimesion.
void
pylith::topology::Field::newSection(const DomainEnum domain,
                                    const int fiberDim,
                                    const int stratum)
{ // newSection
    PYLITH_METHOD_BEGIN;

    // Changing this because cells/vertices are numbered differently in the new scheme
    assert(_dm);
    PetscInt pStart, pEnd;
    PetscErrorCode err;

    switch (domain) {
    case VERTICES_FIELD:
        err = DMPlexGetDepthStratum(_dm, stratum, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        break;
    case CELLS_FIELD:
        err = DMPlexGetHeightStratum(_dm, stratum, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        break;
    case FACES_FIELD:
        err = DMPlexGetHeightStratum(_dm, stratum+1, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        break;
    case POINTS_FIELD:
        err = DMPlexGetChart(_dm, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        break;
    default:
        std::ostringstream msg;
        msg << "Unknown value for DomainEnum: " << domain << "  in Field" << std::endl;
        throw std::logic_error(msg.str());
    }
    newSection(pStart, pEnd, fiberDim);

    PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create PETSc section and set chart and fiber dimesion.
void
pylith::topology::Field::newSection(const PetscInt pStart,
                                    const PetscInt pEnd,
                                    const int fiberDim)
{ // newSection
    PYLITH_METHOD_BEGIN;

    // Changing this because cells/vertices are numbered differently in the new scheme
    assert(_dm);
    PetscSection s = NULL;
    PetscErrorCode err;

    err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(s, pStart, pEnd); PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart; p < pEnd; ++p) {
        err = PetscSectionSetDof(s, p, fiberDim); PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create section given chart.
void
pylith::topology::Field::newSection(const Field& src,
                                    const int fiberDim)
{ // newSection
    PYLITH_METHOD_BEGIN;

    // Clear memory
    clear();
    assert(_dm); assert(src._dm);

    if (fiberDim < 0) {
        std::ostringstream msg;
        msg << "Fiber dimension (" << fiberDim << ") for field '" << _label
            << "' must be nonnegative.";
        throw std::runtime_error(msg.str());
    } // if

    PetscSection srcs = NULL, s = NULL;
    PetscInt pStart, pEnd;
    PetscErrorCode err;

    err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err);

    err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err);
    err = DMGetDefaultSection(src._dm, &srcs); PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetChart(srcs, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(s, pStart, pEnd); PYLITH_CHECK_ERROR(err);
    for (PetscInt p = pStart; p < pEnd; ++p) {
        err = PetscSectionSetDof(s, p, fiberDim); PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // newSection

// ----------------------------------------------------------------------
// Create section with same layout as another section.
void
pylith::topology::Field::cloneSection(const Field& src)
{ // cloneSection
    PYLITH_METHOD_BEGIN;

    const std::string& origLabel = _label;

    // Clear memory
    clear();

    label(origLabel.c_str());

    PetscSection section = src.localSection();
    PetscSection newSection = NULL;
    PetscErrorCode err;

    assert(_dm);
    err = PetscSectionClone(section, &newSection); PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultSection(_dm, newSection); PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&newSection); PYLITH_CHECK_ERROR(err);

    assert(!_localVec);
    err = DMCreateLocalVector(_dm, &_localVec); PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) _localVec,  _label.c_str()); PYLITH_CHECK_ERROR(err);

    // Reuse scatters in clone
    _scatters.clear();
    const scatter_map_type::const_iterator scattersEnd = src._scatters.end();
    for (scatter_map_type::const_iterator s_iter = src._scatters.begin(); s_iter != scattersEnd; ++s_iter) {
        ScatterInfo& sinfo = _scatters[s_iter->first];
        sinfo.dm = NULL;
        sinfo.vector = NULL;

        // Copy DM
        sinfo.dm = s_iter->second.dm;
        err = PetscObjectReference((PetscObject) sinfo.dm); PYLITH_CHECK_ERROR(err);

        // Create global vector.
        err = DMCreateGlobalVector(sinfo.dm, &sinfo.vector); PYLITH_CHECK_ERROR(err);
        err = PetscObjectSetName((PetscObject)sinfo.vector, _label.c_str()); PYLITH_CHECK_ERROR(err);
    } // for

    // Reuse subfields in clone
    _subfields.clear();
    const subfields_type::const_iterator subfieldsEnd = src._subfields.end();
    PetscDS prob = NULL;
    err = DMGetDS(_dm, &prob); PYLITH_CHECK_ERROR(err); assert(prob);
    for (subfields_type::const_iterator s_iter = src._subfields.begin(); s_iter != subfieldsEnd; ++s_iter) {
        SubfieldInfo& sinfo = _subfields[s_iter->first];
        sinfo.description = s_iter->second.description;
        sinfo.index = s_iter->second.index;

        sinfo.fe = s_iter->second.fe;
        PetscFE fe = FieldOps::createFE(sinfo.fe, _dm, _mesh.isSimplex(), sinfo.description.numComponents); assert(fe);
        err = PetscObjectSetName((PetscObject) fe, s_iter->first.c_str()); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetDiscretization(prob, sinfo.index, (PetscObject) fe); PYLITH_CHECK_ERROR(err);
        err = PetscFEDestroy(&fe); PYLITH_CHECK_ERROR(err);

        sinfo.dm = s_iter->second.dm;
        if (sinfo.dm) {
            err = PetscObjectReference((PetscObject) sinfo.dm); PYLITH_CHECK_ERROR(err);
        } // if
    } // for

    PYLITH_METHOD_END;
} // cloneSection

// ----------------------------------------------------------------------
// Clear variables associated with section.
void
pylith::topology::Field::clear(void)
{ // clear
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;

    const scatter_map_type::const_iterator scattersEnd = _scatters.end();
    for (scatter_map_type::iterator s_iter = _scatters.begin(); s_iter != scattersEnd; ++s_iter) {
        err = DMDestroy(&s_iter->second.dm); PYLITH_CHECK_ERROR(err);
        err = VecDestroy(&s_iter->second.vector); PYLITH_CHECK_ERROR(err);
    } // for
    _scatters.clear();

    err = VecDestroy(&_localVec); PYLITH_CHECK_ERROR(err);

    const subfields_type::const_iterator subfieldsEnd = _subfields.end();
    for (subfields_type::iterator s_iter = _subfields.begin(); s_iter != subfieldsEnd; ++s_iter) {
        err = DMDestroy(&s_iter->second.dm); PYLITH_CHECK_ERROR(err);
    } // for

    _dimsOkay = false;

    PYLITH_METHOD_END;
} // clear

// ----------------------------------------------------------------------
// Allocate PETSc section.
void
pylith::topology::Field::allocate(void)
{ // allocate
    PYLITH_METHOD_BEGIN;

    PetscSection s = NULL;
    PetscErrorCode err;

    err = DMGetDefaultSection(_dm, &s); PYLITH_CHECK_ERROR(err); assert(s); // Creates local section
    err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err); // Creates global section
    err = PetscSectionSetUp(s); PYLITH_CHECK_ERROR(err);

    err = VecDestroy(&_localVec); PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(_dm, &_localVec); PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) _localVec,  _label.c_str()); PYLITH_CHECK_ERROR(err);
    err = VecSet(_localVec, 0.0); PYLITH_CHECK_ERROR(err);

    // Create DM for subfields.
    int fields[1];
    for (subfields_type::iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
        SubfieldInfo& info = s_iter->second;
        fields[0] = info.index;
        err = DMDestroy(&info.dm); PYLITH_CHECK_ERROR(err);
        err = DMCreateSubDM(_dm, 1, fields, NULL, &info.dm); PYLITH_CHECK_ERROR(err);
    } // for


    PYLITH_METHOD_END;
} // allocate

// ----------------------------------------------------------------------
// Zero local vector values (including constrained DOF).
void
pylith::topology::Field::zeroLocal(void)
{ // zeroLocal
    PYLITH_METHOD_BEGIN;

    assert(_localVec);
    PetscErrorCode err = VecSet(_localVec, 0.0); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // zeroLocal

// ----------------------------------------------------------------------
// Copy field values and metadata.
void
pylith::topology::Field::copy(const Field& field)
{ // copy
    PYLITH_METHOD_BEGIN;

    // Check compatibility of sections
    const int srcSize = field.chartSize();
    const int dstSize = chartSize();
    if ((field.spaceDim() != spaceDim()) ||
        ( srcSize != dstSize) ) {
        std::ostringstream msg;

        msg << "Cannot copy values from section '" << _label
            << "' to section '" << _label
            << "'. Sections are incompatible.\n"
            << "  Source section size: " << srcSize << "\n"
            << "  Destination section size: " << dstSize;
        throw std::runtime_error(msg.str());
    } // if
    assert(_localVec && field._localVec);

    PetscErrorCode err = VecCopy(field._localVec, _localVec); PYLITH_CHECK_ERROR(err);

    // Update metadata
    label(field._label.c_str());

    PYLITH_METHOD_END;
} // copy

// ----------------------------------------------------------------------
// Dimensionalize field.
void
pylith::topology::Field::dimensionalize(void) const
{ // dimensionalize
    PYLITH_METHOD_BEGIN;

    if (!_dimsOkay) {
        std::ostringstream msg;
        msg << "Cannot dimensionalize field '" << _label
            << "' because the flag has been set to keep field nondimensional.";
        throw std::runtime_error(msg.str());
    } // if

    assert(_localVec);
    const size_t numSubfields = _subfields.size();
    assert(numSubfields > 0);
    if (1 == numSubfields) {
        const PylithReal scale = _subfields.begin()->second.description.scale;
        PetscErrorCode err = VecScale(_localVec, scale); PYLITH_CHECK_ERROR(err);
    } else {
        // Dimensionalize each subfield independently.
        int_array subNumComponents(numSubfields); assert(subNumComponents.size() == numSubfields);
        scalar_array subScales(numSubfields); assert(subScales.size() == numSubfields);
        for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
            const SubfieldInfo& sinfo = s_iter->second;
            const size_t index = sinfo.index; assert(index < numSubfields);
            subScales[index] = sinfo.description.scale;
            subNumComponents[index] = sinfo.description.numComponents;
        } // for

        assert(_dm);
        PylithInt pStart, pEnd;
        PetscErrorCode err = DMPlexGetChart(_dm,  &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        VecVisitorMesh fieldVisitor(*this);
        PylithScalar* fieldArray = fieldVisitor.localArray();
        for (PylithInt p = pStart; p < pEnd; ++p) {
            const size_t dof = fieldVisitor.sectionDof(p);
            if (dof) {
                const PylithInt off = fieldVisitor.sectionOffset(p);
                for (size_t iField = 0, doff = 0; iField < numSubfields; ++iField) {
                    const PylithScalar scale = subScales[iField];
                    const size_t dim = subNumComponents[iField];
                    for (size_t d = 0; d < dim; ++d, ++doff) {
                        assert(doff < dof);
                        fieldArray[off+doff] *= scale;
                    } // for
                } // for
            } // if
        } // for
    } // if/else

    PYLITH_METHOD_END;
} // dimensionalize

// ----------------------------------------------------------------------
// Print field to standard out.
void
pylith::topology::Field::view(const char* label) const
{ // view
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
            std::cout << ", scale: " << sinfo.description.scale << "\n";
        } // for
    } // if

    std::cout << "  dimensionalize flag: " << _dimsOkay << std::endl;

    if (_dm) {
        PetscSection section = NULL;
        PetscMPIInt numProcs, rank;
        PetscErrorCode err;

        err = DMGetDefaultSection(_dm, &section); PYLITH_CHECK_ERROR(err);
        err = DMView(_dm, PETSC_VIEWER_STDOUT_WORLD); PYLITH_CHECK_ERROR(err);
        err = PetscSectionView(section, PETSC_VIEWER_STDOUT_WORLD); PYLITH_CHECK_ERROR(err);
        err = MPI_Comm_size(PetscObjectComm((PetscObject) _dm), &numProcs); PYLITH_CHECK_ERROR(err);
        err = MPI_Comm_rank(PetscObjectComm((PetscObject) _dm), &rank); PYLITH_CHECK_ERROR(err);
        for (PetscInt p = 0; p < numProcs; ++p) {
            err = PetscPrintf(PetscObjectComm((PetscObject) _dm), "Proc %d local vector\n", p); PYLITH_CHECK_ERROR(err);
            if (p == rank) {err = VecView(_localVec, PETSC_VIEWER_STDOUT_SELF); PYLITH_CHECK_ERROR(err); }
            err = PetscBarrier((PetscObject) _dm); PYLITH_CHECK_ERROR(err);
        }
    }

    PYLITH_METHOD_END;
} // view

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// PETSc section view.
void
pylith::topology::Field::createScatter(const Mesh& mesh,
                                       const char* context)
{ // createScatter
    PYLITH_METHOD_BEGIN;

    assert(context);
    PetscErrorCode err;

    const bool createScatterOk = true;
    ScatterInfo& sinfo = _getScatter(context, createScatterOk);
    if (sinfo.dm) {
        assert(sinfo.vector);
        PYLITH_METHOD_END;
    } // if

    err = DMDestroy(&sinfo.dm); PYLITH_CHECK_ERROR(err);
    sinfo.dm = _dm;
    err = PetscObjectReference((PetscObject) sinfo.dm); PYLITH_CHECK_ERROR(err);

    err = VecDestroy(&sinfo.vector); PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(_dm, &sinfo.vector); PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) sinfo.vector, _label.c_str()); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // createScatter

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// PETSc section view. The PETSc vector does not contain constrained
// DOF. Use createScatterWithBC() to include the constrained DOF in
// the PETSc vector.
void
pylith::topology::Field::createScatterWithBC(const Mesh& mesh,
                                             const char* context)
{ // createScatterWithBC
    PYLITH_METHOD_BEGIN;

    assert(context);
    PetscErrorCode err;

    const bool createScatterOk = true;
    ScatterInfo& sinfo = _getScatter(context, createScatterOk);
    if (sinfo.dm) {
        assert(sinfo.vector);
        PYLITH_METHOD_END;
    } // if

    PetscSection section = NULL, newSection = NULL, gsection = NULL;
    PetscDS prob = NULL;
    PetscSF sf = NULL;

    err = DMDestroy(&sinfo.dm); PYLITH_CHECK_ERROR(err);
    err = DMClone(_dm, &sinfo.dm); PYLITH_CHECK_ERROR(err);
    err = DMGetDS(_dm, &prob); PYLITH_CHECK_ERROR(err);
    err = DMSetDS(sinfo.dm, prob); PYLITH_CHECK_ERROR(err);
    err = DMGetDefaultSection(_dm, &section); PYLITH_CHECK_ERROR(err);
    err = PetscSectionClone(section, &newSection); PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultSection(sinfo.dm, newSection); PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&newSection); PYLITH_CHECK_ERROR(err);
    err = DMGetPointSF(sinfo.dm, &sf); PYLITH_CHECK_ERROR(err);
    err = PetscSectionCreateGlobalSection(section, sf, PETSC_TRUE, PETSC_FALSE, &gsection); PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(sinfo.dm, gsection); PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&gsection); PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&sinfo.vector); PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(sinfo.dm, &sinfo.vector); PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) sinfo.vector, _label.c_str()); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // createScatterWithBC

// ----------------------------------------------------------------------
// Create PETSc vector scatter for field. This is used to transfer
// information from the "global" PETSc vector view to the "local"
// PETSc section view. The PETSc vector includes constrained DOF. Use
// createScatter() if constrained DOF should be omitted from the PETSc
// vector.
void
pylith::topology::Field::createScatterWithBC(const Mesh& mesh,
                                             const std::string& labelName,
                                             PetscInt labelValue,
                                             const char* context)
{ // createScatterWithBC
    PYLITH_METHOD_BEGIN;

    assert(context);
    PetscErrorCode err;

    const bool createScatterOk = true;
    ScatterInfo& sinfo = _getScatter(context, createScatterOk);

    // Only create if scatter and scatterVec do not alreay exist.
    if (sinfo.dm) {
        assert(sinfo.vector);
        PYLITH_METHOD_END;
    } // if

    PetscDM dm = mesh.dmMesh(); assert(dm);
    PetscSection section = NULL, newSection = NULL, gsection = NULL, subSection = NULL;
    PetscDS prob = NULL;
    PetscSF sf = NULL;
    PetscDMLabel subpointMap = NULL, subpointMapF = NULL;
    PetscInt dim, dimF, pStart, pEnd, qStart, qEnd, cEnd, cMax, vEnd, vMax;
    err = DMPlexGetHeightStratum(_dm, 0, NULL, &cEnd); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(_dm, 0, NULL, &vEnd); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHybridBounds(_dm, &cMax, NULL, NULL, &vMax); PYLITH_CHECK_ERROR(err);
    PetscInt excludeRanges[4] = {cMax, cEnd, vMax, vEnd};
    PetscInt numExcludes = (cMax >= 0 ? 1 : 0) + (vMax >= 0 ? 1 : 0);

    err = DMGetDefaultSection(_dm, &section); PYLITH_CHECK_ERROR(err);
    err = DMGetDimension(dm,  &dim); PYLITH_CHECK_ERROR(err);
    err = DMGetDimension(_dm, &dimF); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetChart(dm,  &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetChart(_dm, &qStart, &qEnd); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSubpointMap(dm,  &subpointMap); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSubpointMap(_dm, &subpointMapF); PYLITH_CHECK_ERROR(err);
    if (((dim != dimF) || ((pEnd-pStart) < (qEnd-qStart))) && subpointMap && !subpointMapF) {
        const PetscInt *ind = NULL;
        PetscIS subpointIS = NULL;
        PetscInt n = 0, q = 0;

        err = PetscSectionGetChart(section, &qStart, &qEnd); PYLITH_CHECK_ERROR(err);
        err = DMPlexCreateSubpointIS(dm, &subpointIS); PYLITH_CHECK_ERROR(err);
        if (subpointIS) {
            err = ISGetLocalSize(subpointIS, &n); PYLITH_CHECK_ERROR(err);
            err = ISGetIndices(subpointIS, &ind); PYLITH_CHECK_ERROR(err);
        } // if
        err = PetscSectionCreate(mesh.comm(), &subSection); PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetChart(subSection, pStart, pEnd); PYLITH_CHECK_ERROR(err);
        for (q = qStart; q < qEnd; ++q) {
            PetscInt dof, off, p;

            err = PetscSectionGetDof(section, q, &dof); PYLITH_CHECK_ERROR(err);
            if (dof) {
                err = PetscFindInt(q, n, ind, &p); PYLITH_CHECK_ERROR(err);
                if ((p >= pStart) && (p < pEnd)) {
                    err = PetscSectionSetDof(subSection, p, dof); PYLITH_CHECK_ERROR(err);
                    err = PetscSectionGetOffset(section, q, &off); PYLITH_CHECK_ERROR(err);
                    err = PetscSectionSetOffset(subSection, p, off); PYLITH_CHECK_ERROR(err);
                } // if
            } // if
        } // for
        if (subpointIS) {
            err = ISRestoreIndices(subpointIS, &ind); PYLITH_CHECK_ERROR(err);
            err = ISDestroy(&subpointIS); PYLITH_CHECK_ERROR(err);
        } // if
          /* No need to setup section */
        section = subSection;
        /* There are no excludes for surface meshes */
        numExcludes = 0;
    } // if

    err = DMDestroy(&sinfo.dm); PYLITH_CHECK_ERROR(err);
    err = DMClone(_dm, &sinfo.dm); PYLITH_CHECK_ERROR(err);
    err = DMGetDS(_dm, &prob); PYLITH_CHECK_ERROR(err);
    err = DMSetDS(sinfo.dm, prob); PYLITH_CHECK_ERROR(err);
    err = PetscSectionClone(section, &newSection); PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultSection(sinfo.dm, newSection); PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&newSection); PYLITH_CHECK_ERROR(err);
    err = DMGetPointSF(sinfo.dm, &sf); PYLITH_CHECK_ERROR(err);
    if (labelName.empty()) {
        err = PetscSectionCreateGlobalSectionCensored(section, sf, PETSC_TRUE, numExcludes, excludeRanges, &gsection);
        PYLITH_CHECK_ERROR(err);
    } else {
        DMLabel label;

        err = DMGetLabel(sinfo.dm, labelName.c_str(), &label); PYLITH_CHECK_ERROR(err);
        err = PetscSectionCreateGlobalSectionLabel(section, sf, PETSC_TRUE, label, labelValue, &gsection);
        PYLITH_CHECK_ERROR(err);
    } // if/else
    err = DMSetDefaultGlobalSection(sinfo.dm, gsection); PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&gsection); PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&sinfo.vector); PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(sinfo.dm, &sinfo.vector); PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) sinfo.vector, _label.c_str()); PYLITH_CHECK_ERROR(err);

    err = PetscSectionDestroy(&subSection); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // createScatterWithBC

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
PetscVec
pylith::topology::Field::scatterVector(const char* context)
{ // scatterVector
    PYLITH_METHOD_BEGIN;

    ScatterInfo& sinfo = _getScatter(context);

    PYLITH_METHOD_RETURN(sinfo.vector);
} // scatterVector

// ----------------------------------------------------------------------
// Get PETSc vector associated with field.
const PetscVec
pylith::topology::Field::scatterVector(const char* context) const
{ // scatterVector
    PYLITH_METHOD_BEGIN;

    const ScatterInfo& sinfo = _getScatter(context);

    PYLITH_METHOD_RETURN(sinfo.vector);
} // scatterVector

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
// PETSc vector view of the field.
void
pylith::topology::Field::scatterLocalToContext(const char* context,
                                               InsertMode mode) const
{ // scatterLocalToContext
    PYLITH_METHOD_BEGIN;

    assert(context);

    const ScatterInfo& sinfo = _getScatter(context);
    PetscErrorCode err;
    if (sinfo.dm) {
        err = DMLocalToGlobalBegin(sinfo.dm, _localVec, mode, sinfo.vector); PYLITH_CHECK_ERROR(err);
        err = DMLocalToGlobalEnd(sinfo.dm, _localVec, mode, sinfo.vector); PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // scatterLocalToContext

// ----------------------------------------------------------------------
// Scatter section information across processors to update the
// PETSc vector view of the field.
void
pylith::topology::Field::scatterLocalToVector(const PetscVec vector,
                                              InsertMode mode) const
{ // scatterLocalToVector
    PYLITH_METHOD_BEGIN;

    assert(vector);

    PetscErrorCode err;
    if (_dm) {
        err = DMLocalToGlobalBegin(_dm, _localVec, mode, vector); PYLITH_CHECK_ERROR(err);
        err = DMLocalToGlobalEnd(_dm, _localVec, mode, vector); PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // scatterLocalToVector

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
void
pylith::topology::Field::scatterContextToLocal(const char* context,
                                               InsertMode mode) const
{ // scatterContextToLocal
    PYLITH_METHOD_BEGIN;

    assert(context);

    const ScatterInfo& sinfo = _getScatter(context);
    PetscErrorCode err;
    if (sinfo.dm) {
        err = DMGlobalToLocalBegin(sinfo.dm, sinfo.vector, mode, _localVec); PYLITH_CHECK_ERROR(err);
        err = DMGlobalToLocalEnd(sinfo.dm, sinfo.vector, mode, _localVec); PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // scatterContextToLocal

// ----------------------------------------------------------------------
// Scatter PETSc vector information across processors to update the
// section view of the field.
void
pylith::topology::Field::scatterVectorToLocal(const PetscVec vector,
                                              InsertMode mode) const
{ // scatterVectorToLocal
    PYLITH_METHOD_BEGIN;

    assert(vector);

    PetscErrorCode err;
    if (_dm) {
        err = DMGlobalToLocalBegin(_dm, vector, mode, _localVec); PYLITH_CHECK_ERROR(err);
        err = DMGlobalToLocalEnd(_dm, vector, mode, _localVec); PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // scatterVectorToLocal

// ----------------------------------------------------------------------
// Get scatter for given context.
pylith::topology::Field::ScatterInfo&
pylith::topology::Field::_getScatter(const char* context,
                                     const bool createOk)
{ // _getScatter
    PYLITH_METHOD_BEGIN;

    assert(context);

    bool isNewScatter = _scatters.find(context) == _scatters.end();

    // Synchronize creation of scatter (empty sections may have
    // leftover, reusable scatters that need to be cleared out).
    int numNewScatterLocal = (isNewScatter) ? 1 : 0;
    int numNewScatter = 0;
    MPI_Allreduce(&numNewScatterLocal, &numNewScatter, 1, MPI_INT, MPI_MAX, _mesh.comm());
    if (numNewScatter && !isNewScatter) {
        // remove old scatter
        ScatterInfo& sinfo = _scatters[context];
        PetscErrorCode err;
        err = DMDestroy(&sinfo.dm); PYLITH_CHECK_ERROR(err);
        err = VecDestroy(&sinfo.vector); PYLITH_CHECK_ERROR(err);

        _scatters.erase(context);
        isNewScatter = true;
    } // if

    if (isNewScatter && !createOk) {
        std::ostringstream msg;
        msg << "Scatter for context '" << context << "' does not exist for field '" << label() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    ScatterInfo& sinfo = _scatters[context];
    if (isNewScatter) {
        sinfo.dm = NULL;
        sinfo.vector = NULL;
    } // if
    assert(_scatters.find(context) != _scatters.end());

    PYLITH_METHOD_RETURN(sinfo);
} // _getScatter

// ----------------------------------------------------------------------
// Get scatter for given context.
const pylith::topology::Field::ScatterInfo&
pylith::topology::Field::_getScatter(const char* context) const
{ // _getScatter
    PYLITH_METHOD_BEGIN;

    assert(context);

    const scatter_map_type::const_iterator s_iter =
        _scatters.find(context);
    if (s_iter == _scatters.end()) {
        std::ostringstream msg;
        msg << "Scatter for context '" << context << "' does not exist for field '" << label() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(s_iter->second);
} // _getScatter

// ----------------------------------------------------------------------
// Add subfield to current field (inteface to use from SWIG).
void
pylith::topology::Field::subfieldAdd(const char *name,
                                     const VectorFieldEnum fieldType,
                                     const char* components[],
                                     const int numComponents,
                                     const double scale,
                                     const int basisOrder,
                                     const int quadOrder,
                                     const bool isBasisContinuous,
                                     const SpaceEnum feSpace)
{ // subfieldAdd
    assert(numComponents > 0);

    Description description;
    description.label = name;
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
    discretization.isBasisContinuous = isBasisContinuous;
    discretization.feSpace = feSpace;

    this->subfieldAdd(description, discretization);
} // subfieldAdd

// ----------------------------------------------------------------------
// Add subfield.
void
pylith::topology::Field::subfieldAdd(const Description& description,
                                     const Discretization& discretization)
{ // subfieldAdd
    PYLITH_METHOD_BEGIN;

    assert(0 == _subfields.count(description.label));

    // Keep track of name/components until setup
    SubfieldInfo info;
    info.description = description;
    info.fe = discretization;
    info.index = _subfields.size(); // Indices match order added.
    info.dm = NULL;
    _subfields[description.label] = info;

    PYLITH_METHOD_END;
} // subfieldAdd

// ----------------------------------------------------------------------
void
pylith::topology::Field::subfieldsSetup(void)
{ // subfieldsSetup
    PYLITH_METHOD_BEGIN;

    assert(_dm);

    // Setup section now that we know the total number of sub-fields and components.
    PetscDS prob = NULL;
    PetscErrorCode err;

    err = DMGetDS(_dm, &prob); PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultSection(_dm, NULL); PYLITH_CHECK_ERROR(err); // :TODO: @brad Remove when using PetscDS for all fields.
    err = DMSetNumFields(_dm, _subfields.size()); PYLITH_CHECK_ERROR(err);

    for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
        const char* sname = s_iter->first.c_str();
        const SubfieldInfo& sinfo = s_iter->second;

        PetscFE fe = FieldOps::createFE(sinfo.fe, _dm, _mesh.isSimplex(), sinfo.description.numComponents); assert(fe);
        err = PetscObjectSetName((PetscObject) fe, sname); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetDiscretization(prob, sinfo.index, (PetscObject) fe); PYLITH_CHECK_ERROR(err);
        err = PetscFEDestroy(&fe); PYLITH_CHECK_ERROR(err);
    } // for
    err = PetscDSSetUp(prob); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // subfieldsSetup

// ----------------------------------------------------------------------
void
pylith::topology::Field::subfieldSetDof(const char *name,
                                        const DomainEnum domain,
                                        int fiberDim)
{ // subfieldSetDof
    PYLITH_METHOD_BEGIN;

    // :TEMPORARY: Obsolete when using PetscDS

    PetscInt pStart, pEnd;
    PetscErrorCode err;

    assert(_dm);
    switch (domain) {
    case VERTICES_FIELD:
        err = DMPlexGetDepthStratum(_dm, 0, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        break;
    case CELLS_FIELD:
        err = DMPlexGetHeightStratum(_dm, 0, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        break;
    case FACES_FIELD:
        err = DMPlexGetHeightStratum(_dm, 1, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        break;
    case POINTS_FIELD:
        err = DMPlexGetChart(_dm, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        break;
    default:
        std::ostringstream msg;
        msg << "Unknown value for DomainEnum: " << domain << "  in Field" << std::endl;
        throw std::logic_error(msg.str());
    } // switch
    PetscSection section = NULL;
    err = DMGetDefaultSection(_dm, &section); PYLITH_CHECK_ERROR(err); assert(section);
    const int iField = _subfields[name].index;
    for (PetscInt p = pStart; p < pEnd; ++p) {
        //err = PetscSectionAddDof(section, p, fiberDim);PYLITH_CHECK_ERROR(err); // Future use
        err = PetscSectionSetFieldDof(section, p, iField, fiberDim); PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // subfieldSetDof

// ----------------------------------------------------------------------
// Does field have given auxiliary subfield?
bool
pylith::topology::Field::hasSubfield(const char* name) const
{ // hasSubfield
    PYLITH_METHOD_BEGIN;

    subfields_type::const_iterator iter = _subfields.find(name);
    PYLITH_METHOD_RETURN(_subfields.end() != iter);
} // hasSubfield


// ----------------------------------------------------------------------
// Get names of subfields.
pylith::string_vector
pylith::topology::Field::subfieldNames(void) const
{  // subfieldNames
    PYLITH_METHOD_BEGIN;

    const size_t numSubfields = _subfields.size();
    pylith::string_vector names(numSubfields);
    for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
        const SubfieldInfo& sinfo = s_iter->second;
        names[sinfo.index] = s_iter->first;
    } // for

    PYLITH_METHOD_RETURN(pylith::string_vector(names));
}  // subfieldNames


// ----------------------------------------------------------------------
// Get metadata for subfield.
const pylith::topology::Field::SubfieldInfo&
pylith::topology::Field::subfieldInfo(const char* name) const
{ // subfieldInfo
    PYLITH_METHOD_BEGIN;

    subfields_type::const_iterator iter = _subfields.find(name);
    if (_subfields.end() == iter) {
        std::ostringstream msg;
        msg << "Could not find subfield '" << name << "' in field '" << label() << "'.\n"
            << "Available subfields:";
        for (subfields_type::const_iterator s_iter = _subfields.begin(); s_iter != _subfields.end(); ++s_iter) {
            msg << " '" << s_iter->first << "'";
        } // for
        msg << std::endl;

        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_RETURN(iter->second);
} // subfieldInfo


// ----------------------------------------------------------------------
// Copy subfield values to field.
void
pylith::topology::Field::copySubfield(const Field& field,
                                      const char* name)
{ // copySubfield
    PYLITH_METHOD_BEGIN;

    // Check compatibility of sections
    PYLITH_JOURNAL_ERROR(":TODO: @brad This is an insufficient test of field compatibility.");
    const int srcSize = field.chartSize();
    const int dstSize = chartSize();
    if (dstSize != srcSize) {
        _extractSubfield(field, name);
    } // if
    assert(_localVec && field._localVec);

    const SubfieldInfo& subfieldInfo = const_cast<Field&>(field)._subfields[name];
    const int subfieldIndex = subfieldInfo.index; assert(subfieldIndex >= 0);

    label(subfieldInfo.description.label.c_str()); // Use method to insure propagation to subsidiary objects

    PetscErrorCode err;
    const PetscSection& fieldSection = field.localSection();
    const PetscSection& subfieldSection = this->localSection();

    PetscInt pStart, pEnd;
    err = PetscSectionGetChart(subfieldSection, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);

    // Copy values from field
    PylithScalar* subfieldArray = NULL;
    PylithScalar* fieldArray = NULL;
    err = VecGetArray(this->_localVec, &subfieldArray); PYLITH_CHECK_ERROR(err);
    err = VecGetArray(field._localVec, &fieldArray); PYLITH_CHECK_ERROR(err);
    for (PetscInt p = pStart; p < pEnd; ++p) {
        PetscInt fdof, foff, sdof, soff;

        err = PetscSectionGetDof(subfieldSection, p, &sdof); PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetOffset(subfieldSection, p, &soff); PYLITH_CHECK_ERROR(err);

        err = PetscSectionGetFieldDof(fieldSection, p, subfieldIndex, &fdof); PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetFieldOffset(fieldSection, p, subfieldIndex, &foff); PYLITH_CHECK_ERROR(err);

        assert(fdof == sdof);
        for (PetscInt d = 0; d < fdof; ++d) {
            subfieldArray[soff+d] = fieldArray[foff+d];
        } // for
    } // for
    err = VecRestoreArray(field._localVec, &fieldArray); PYLITH_CHECK_ERROR(err);
    err = VecRestoreArray(this->_localVec, &subfieldArray); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // copySubfield


// ----------------------------------------------------------------------
// Extract subfield.
void
pylith::topology::Field::_extractSubfield(const Field& field,
                                          const char* name)
{ // _extractSubfield
    PYLITH_METHOD_BEGIN;

    clear();

    const SubfieldInfo& subfieldInfo = const_cast<Field&>(field)._subfields[name];
    const int subfieldIndex = subfieldInfo.index; assert(subfieldIndex >= 0);

    PetscErrorCode err;
    PetscIS subfieldIS = NULL;
    const int numSubfields = 1;
    int indicesSubfield[1];
    indicesSubfield[0] = subfieldIndex;
    err = DMDestroy(&_dm); PYLITH_CHECK_ERROR(err);
    if (subfieldInfo.dm) {
        PetscSection s;
        err = DMClone(subfieldInfo.dm, &_dm); PYLITH_CHECK_ERROR(err); assert(_dm);
        err = DMGetDefaultSection(subfieldInfo.dm, &s); PYLITH_CHECK_ERROR(err);
        err = DMSetDefaultSection(_dm, s); PYLITH_CHECK_ERROR(err);
    } else {
        err = DMCreateSubDM(field.dmMesh(), numSubfields, indicesSubfield, &subfieldIS, &_dm); PYLITH_CHECK_ERROR(err);
        assert(_dm);
    } // if/else
    err = ISDestroy(&subfieldIS); PYLITH_CHECK_ERROR(err);

    this->subfieldAdd(subfieldInfo.description, subfieldInfo.fe);
    this->subfieldsSetup();

    err = DMCreateLocalVector(_dm, &_localVec); PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) _localVec,  _label.c_str()); PYLITH_CHECK_ERROR(err);

    // Setup section
    const PetscSection& fieldSection = field.localSection();
    PetscSection subfieldSection = NULL;
    err = DMGetDefaultSection(_dm, &subfieldSection); PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(_dm, NULL); PYLITH_CHECK_ERROR(err);

    PetscInt pStart = -1, pEnd = -1;
    err = PetscSectionGetChart(fieldSection, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
    err = PetscSectionSetChart(subfieldSection, pStart, pEnd); PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart; p < pEnd; ++p) {
        PetscInt dof;
        err = PetscSectionGetFieldDof(fieldSection, p, subfieldIndex, &dof); PYLITH_CHECK_ERROR(err);
        if (dof > 0) {
            err = PetscSectionSetDof(subfieldSection, p, dof); PYLITH_CHECK_ERROR(err);
            err = PetscSectionSetFieldDof(subfieldSection, p, 0, dof); PYLITH_CHECK_ERROR(err);

            err = PetscSectionGetFieldConstraintDof(fieldSection, p, subfieldIndex, &dof); PYLITH_CHECK_ERROR(err);
            err = PetscSectionSetConstraintDof(subfieldSection, p, dof); PYLITH_CHECK_ERROR(err);
        } // if
    } // for
    allocate();

    for (PetscInt p = pStart; p < pEnd; ++p) {
        PetscInt dof;
        const PetscInt* indices = NULL;
        err = PetscSectionGetConstraintDof(subfieldSection, p, &dof); PYLITH_CHECK_ERROR(err);
        if (dof > 0) {
            err = PetscSectionGetFieldConstraintIndices(fieldSection, p, subfieldIndex, &indices); PYLITH_CHECK_ERROR(
                err);
            err = PetscSectionSetConstraintIndices(subfieldSection, p, indices); PYLITH_CHECK_ERROR(err);
        } // if
    } // for

    PYLITH_METHOD_END;
} // _extractSubField

// End of file
