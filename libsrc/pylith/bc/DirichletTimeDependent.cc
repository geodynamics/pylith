// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "DirichletTimeDependent.hh" // implementation of object methods

#include "DirichletAuxiliaryFactory.hh" // USES DirichletAuxiliaryFactory

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

extern "C" {
    #include "pylith/fekernels/timedependentbc.h"
}

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
const char* pylith::bc::DirichletTimeDependent::_pyreComponent = "dirichlettimedependent";

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletTimeDependent::DirichletTimeDependent(void) :
    _dbTimeHistory(NULL),
    _useInitial(true),
    _useRate(false),
    _useTimeHistory(false)
{ // constructor
    PyreComponent::name(_pyreComponent);
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletTimeDependent::~DirichletTimeDependent(void)
{ // destructor
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletTimeDependent::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    DirichletNew::deallocate();
    _dbTimeHistory = NULL; // :KLUDGE: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set time history database.
void
pylith::bc::DirichletTimeDependent::dbTimeHistory(spatialdata::spatialdb::TimeHistory* th)
{ // dbTimeHistory
    _dbTimeHistory = th;
} // dbTimeHistory


// ----------------------------------------------------------------------
// Get time history database.
const spatialdata::spatialdb::TimeHistory*
pylith::bc::DirichletTimeDependent::dbTimeHistory(void)
{ // dbTimeHistory
    return _dbTimeHistory;
} // dbTimeHistory

// ----------------------------------------------------------------------
// Use initial value term in time history expression.
void
pylith::bc::DirichletTimeDependent::useInitial(const bool value)
{ // useInitial
    PYLITH_COMPONENT_DEBUG("useInitial(value="<<value<<")");

    _useInitial = value;
} // useInitial


// ----------------------------------------------------------------------
// Get flag associated with using initial value term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useInitial(void) const
{ // useInitial
    return _useInitial;
} // useInitial


// ----------------------------------------------------------------------
// Use rate value term in time history expression.
void
pylith::bc::DirichletTimeDependent::useRate(const bool value)
{ // useRate
    PYLITH_COMPONENT_DEBUG("useRate(value="<<value<<")");

    _useRate = value;
} // useRate


// ----------------------------------------------------------------------
// Get flag associated with using rate value term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useRate(void) const
{ // useRate
    return _useRate;
} // useRate


// ----------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::bc::DirichletTimeDependent::useTimeHistory(const bool value)
{ // useTimeHistory
    PYLITH_COMPONENT_DEBUG("useTimeHistory(value="<<value<<")");

    _useTimeHistory = value;
} // useTimeHistory


// ----------------------------------------------------------------------
// Get flag associated with using time history term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useTimeHistory(void) const
{ // useTimeHistory
    return _useTimeHistory;
} // useTimeHistory


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::bc::DirichletTimeDependent::prestep(const double t,
                                            const double dt)
{ // prestep
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_prestep(t="<<t<<", dt="<<dt<<")");

    if (_useTimeHistory) {
        assert(_normalizer);
        assert(_auxFields);

        const PylithScalar timeScale = _normalizer->timeScale();

        PetscErrorCode err;

        PetscSection auxFieldsSection = _auxFields->localSection(); assert(auxFieldsSection);
        PetscInt pStart = 0, pEnd = 0;
        err = PetscSectionGetChart(auxFieldsSection, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        pylith::topology::VecVisitorMesh auxFieldsVisitor(*_auxFields);
        PetscScalar* auxFieldsArray = auxFieldsVisitor.localArray(); assert(auxFieldsArray);

        const PylithInt numComponents = _description.numComponents;
        PetscInt offTH = 0;
        if (_useInitial) {offTH += numComponents;}
        if (_useRate) {offTH += numComponents + 1;}
        const PetscInt offStartTime = offTH + numComponents;
        const PetscInt offValue = offStartTime + 1;

        // Loop over all points in section.
        for (PetscInt p = pStart; p < pEnd; ++p) {
            // Skip points without values in section.
            if (!auxFieldsVisitor.sectionDof(p)) {continue;}

            // Get offset for point.
            const PetscInt off = auxFieldsVisitor.sectionOffset(p);

            // Get starting time and compute relative time for point.
            const PylithScalar tStart = auxFieldsArray[off+offStartTime];
            const PylithScalar tRel = t - tStart;

            // Query time history for value (normalized amplitude).
            PylithScalar value = 0.0;
            if (tRel >= 0.0) {
                PylithScalar tDim = tRel * timeScale;
                const int err = _dbTimeHistory->query(&value, tDim);
                if (err) {
                    std::ostringstream msg;
                    msg << "Error querying for time '" << tDim << "' in time history database '" << _dbTimeHistory->label() << "'.";
                    throw std::runtime_error(msg.str());
                } // if
            } // if
              // Update value (normalized amplitude) in auxiliary field.
            auxFieldsArray[off+offValue] = value;
        } // for

    } // if

    PYLITH_METHOD_END;
} // prestep


// ----------------------------------------------------------------------
// Setup auxiliary subfields (discretization and query fns).
void
pylith::bc::DirichletTimeDependent::_auxFieldsSetup(const pylith::topology::Field& solution)
{ // _auxFieldsSetup
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldsSetup(solution="<<solution.label()<<")");

    DirichletAuxiliaryFactory factory(*this, solution, _normalizer->timeScale());

    // :ATTENTION: The order for subfieldAdd() must match the order of the auxiliary fields in the FE kernels.

    if (_useInitial) {
        factory.initialAmplitude();
    } // if
    if (_useRate) {
        factory.rateAmplitude();
        factory.rateStartTime();
    } // _useRate
    if (_useTimeHistory) {
        factory.timeHistoryAmplitude();
        factory.timeHistoryStartTime();
        factory.timeHistoryValue();
    } // _useTimeHistory

    PYLITH_METHOD_END;
}     // _auxFieldsSetup


// ----------------------------------------------------------------------
// Set kernels for setting Dirhclet values.
void
pylith::bc::DirichletTimeDependent::_setFEKernelsConstraint(const pylith::topology::Field& solution)
{ // _setFEKernelsConstraint
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsConstraint(solution="<<solution.label()<<")");

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    const bool isScalarField = _description.vectorFieldType == pylith::topology::Field::SCALAR;

    const int bitInitial = _useInitial ? 0x1 : 0x0;
    const int bitRate = _useRate ? 0x2 : 0x0;
    const int bitTimeHistory = _useTimeHistory ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;
    switch (bitUse) {
    case 0x1:
        _bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_initial_scalar : pylith_fekernels_TimeDependentBC_initial_vector;
        break;
    case 0x2:
        _bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_rate_scalar : pylith_fekernels_TimeDependentBC_rate_vector;
        break;
    case 0x4:
        _bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_timeHistory_scalar : pylith_fekernels_TimeDependentBC_timeHistory_vector;
        break;
    case 0x3:
        _bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_initialRate_scalar : pylith_fekernels_TimeDependentBC_initialRate_vector;
        break;
    case 0x5:
        _bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_initialTimeHistory_scalar : pylith_fekernels_TimeDependentBC_initialTimeHistory_vector;
        break;
    case 0x6:
        _bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_rateTimeHistory_scalar : pylith_fekernels_TimeDependentBC_rateTimeHistory_vector;
        break;
    case 0x7:
        _bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_initialRateTimeHistory_scalar : pylith_fekernels_TimeDependentBC_initialRateTimeHistory_vector;
        break;
    case 0x0:
        PYLITH_COMPONENT_WARNING("Dirichlet BC provides no constraints.");
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown combination of flags for Dirichlet BC terms (useInitial="<<_useInitial<<", useRate="<<_useRate<<", useTimeHistory="<<_useTimeHistory<<").");
        throw std::logic_error("Unknown combination of flags for Dirichlet BC terms.");
    } // switch


    PYLITH_METHOD_END;
} // _setFEKernelsConstraint

// End of file
