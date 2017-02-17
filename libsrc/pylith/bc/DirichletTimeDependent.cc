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

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

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
    PYLITH_JOURNAL_DEBUG("useInitial(value="<<value<<")");

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
    PYLITH_JOURNAL_DEBUG("useRate(value="<<value<<")");

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
    PYLITH_JOURNAL_DEBUG("useTimeHistory(value="<<value<<")");

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
    PYLITH_JOURNAL_DEBUG("_prestep(t="<<t<<", dt="<<dt<<")");

    if (_useTimeHistory) {
        assert(_normalizer);
        assert(_auxFields);

        const PylithScalar timeScale = _normalizer->timeScale();

        PetscErrorCode err;

        PetscSection auxFieldsSection = _auxFields->localSection(); assert(auxFieldsSection);
        PetscInt pStart=0, pEnd=0;
        err = PetscSectionGetChart(auxFieldsSection, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        pylith::topology::VecVisitorMesh auxFieldsVisitor(*_auxFields);
        PetscScalar* auxFieldsArray = auxFieldsVisitor.localArray(); assert(auxFieldsArray);

        // Compute offset of time history subfields in auxiliary field.
        // :ASSUMPTION: Constrained field is a scalar or vector field.
        const PetscInt numComponents = (_vectorFieldType == pylith::topology::Field::VECTOR) ? _spaceDim : 1;
        PetscInt offTH = 0;
        if (_useInitial) offTH += numComponents;
        if (_useRate) offTH += numComponents + 1;
        const PetscInt offStartTime = offTH + numComponents;
        const PetscInt offValue = offStartTime + 1;

        // Loop over all points in section.
        for (PetscInt p=pStart; p < pEnd; ++p) {
            // Skip points without values in section.
            if (!auxFieldsVisitor.sectionDof(p)) continue;

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
pylith::bc::DirichletTimeDependent::_auxFieldsSetup(void)
{ // _auxFieldsSetup
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_auxFieldsSetup()");

    // Set subfields in auxiliary fields.
    assert(_normalizer);
    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal timeScale = _normalizer->timeScale();
    const PylithReal velocityScale = lengthScale / timeScale;

    // :ASSUMPTION: Constrained field is a scalar or vector field.
    const bool isVector = _vectorFieldType == pylith::topology::Field::VECTOR;
    const int numComponents = (isVector) ? _spaceDim : 1;

    // :ATTENTION: The order for subfieldAdd() must match the order of the auxiliary fields in the FE kernels.

    // Initial amplitude
    if (_useInitial) {
        const pylith::topology::Field::DiscretizeInfo& initialAmplitudeFEInfo = this->auxFieldDiscretization("initial_amplitude");
        if (isVector) {
            const char* componentNames[3] = {"initial_amplitude_x", "initial_amplitude_y", "initial_amplitude_z"};
            _auxFields->subfieldAdd("initial_amplitude", componentNames, numComponents, _vectorFieldType, initialAmplitudeFEInfo.basisOrder, initialAmplitudeFEInfo.quadOrder, initialAmplitudeFEInfo.isBasisContinuous, lengthScale);
        } else {
            const char* componentNames[1] = {"initial_amplitude"};
            _auxFields->subfieldAdd("initial_amplitude", componentNames, numComponents, _vectorFieldType, initialAmplitudeFEInfo.basisOrder, initialAmplitudeFEInfo.quadOrder, initialAmplitudeFEInfo.isBasisContinuous, lengthScale);
        } // if/else
        _auxFieldsQuery->queryFn("initial_amplitude", pylith::topology::FieldQuery::dbQueryGeneric);
    } // if


    // Rate amplitude and start time.
    if (_useRate) {
        const pylith::topology::Field::DiscretizeInfo& rateAmplitudeFEInfo = this->auxFieldDiscretization("rate_amplitude");
        if (isVector) {
            const char* componentNames[3] = {"rate_amplitude_x", "rate_amplitude_y", "rate_amplitude_z"};
            _auxFields->subfieldAdd("rate_amplitude", componentNames, numComponents, _vectorFieldType, rateAmplitudeFEInfo.basisOrder, rateAmplitudeFEInfo.quadOrder, rateAmplitudeFEInfo.isBasisContinuous, velocityScale);
        } else {
            const char* componentNames[1] = {"rate_amplitude"};
            _auxFields->subfieldAdd("rate_amplitude", componentNames, numComponents, _vectorFieldType, rateAmplitudeFEInfo.basisOrder, rateAmplitudeFEInfo.quadOrder, rateAmplitudeFEInfo.isBasisContinuous, velocityScale);
        } // if/else
        _auxFieldsQuery->queryFn("rate_amplitude", pylith::topology::FieldQuery::dbQueryGeneric);

        const char* startNames[1] = {"rate_start"};
        const pylith::topology::Field::DiscretizeInfo& rateStartFEInfo = this->auxFieldDiscretization("rate_start");
        _auxFields->subfieldAdd("rate_start", startNames, 1, pylith::topology::Field::SCALAR, rateStartFEInfo.basisOrder, rateStartFEInfo.quadOrder, rateStartFEInfo.isBasisContinuous, timeScale);
        _auxFieldsQuery->queryFn("rate_start", pylith::topology::FieldQuery::dbQueryGeneric);
    } // if


    // Time history amplitude and start time.
    if (_useTimeHistory) {
        const pylith::topology::Field::DiscretizeInfo& timeHistoryAmplitudeFEInfo = this->auxFieldDiscretization("time_history_amplitude");
        if (isVector) {
            const char* componentNames[3] = {"time_history_amplitude_x", "time_history_amplitude_y", "time_history_amplitude_z"};
            _auxFields->subfieldAdd("time_history_amplitude", componentNames, numComponents, _vectorFieldType, timeHistoryAmplitudeFEInfo.basisOrder, timeHistoryAmplitudeFEInfo.quadOrder, timeHistoryAmplitudeFEInfo.isBasisContinuous, lengthScale);
        } else {
            const char* componentNames[1] = {"time_history_amplitude"};
            _auxFields->subfieldAdd("time_history_amplitude", componentNames, numComponents, _vectorFieldType, timeHistoryAmplitudeFEInfo.basisOrder, timeHistoryAmplitudeFEInfo.quadOrder, timeHistoryAmplitudeFEInfo.isBasisContinuous, lengthScale);
        } // if/else
        _auxFieldsQuery->queryFn("time_history_amplitude", pylith::topology::FieldQuery::dbQueryGeneric);

        const char* startNames[1] = {"time_history_start"};
        const pylith::topology::Field::DiscretizeInfo& timeHistoryStartFEInfo = this->auxFieldDiscretization("time_history_start");
        _auxFields->subfieldAdd("time_history_start", startNames, 1, pylith::topology::Field::SCALAR, timeHistoryStartFEInfo.basisOrder, timeHistoryStartFEInfo.quadOrder, timeHistoryStartFEInfo.isBasisContinuous, timeScale);
        _auxFieldsQuery->queryFn("time_history_start", pylith::topology::FieldQuery::dbQueryGeneric);

        // Field to hold value from query of time history database.
        const char* valueNames[1] = {"time_history_value"};
        _auxFields->subfieldAdd("time_history_value", valueNames, 1, pylith::topology::Field::SCALAR, timeHistoryAmplitudeFEInfo.basisOrder, timeHistoryAmplitudeFEInfo.quadOrder, timeHistoryAmplitudeFEInfo.isBasisContinuous, 1.0);
        _auxFieldsQuery->queryFn("time_history_value", NULL);
    } // if


    PYLITH_METHOD_END;
}     // _auxFieldsSetup


// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::bc::DirichletTimeDependent::_setFEKernelsConstraint(const topology::Field& solution)
{ // _setFEKernelsConstraint
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setFEKernelsConstraint(solution="<<solution.label()<<")");

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    const bool isScalarField = _vectorFieldType == pylith::topology::Field::SCALAR;

    PetscPointFunc bcKernel = NULL;
    const int bitInitial = _useInitial ? 0x1 : 0x0;
    const int bitRate = _useRate ? 0x2 : 0x0;
    const int bitTimeHistory = _useTimeHistory ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;
    switch (bitUse) {
    case 0x1:
        bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_initial_scalar : pylith_fekernels_TimeDependentBC_initial_vector;
        break;
    case 0x2:
        bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_rate_scalar : pylith_fekernels_TimeDependentBC_rate_vector;
        break;
    case 0x4:
        bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_timeHistory_scalar : pylith_fekernels_TimeDependentBC_timeHistory_vector;
        break;
    case 0x3:
        bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_initialRate_scalar : pylith_fekernels_TimeDependentBC_initialRate_vector;
        break;
    case 0x5:
        bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_initialTimeHistory_scalar : pylith_fekernels_TimeDependentBC_initialTimeHistory_vector;
        break;
    case 0x6:
        bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_rateTimeHistory_scalar : pylith_fekernels_TimeDependentBC_rateTimeHistory_vector;
        break;
    case 0x7:
        bcKernel = (isScalarField) ? pylith_fekernels_TimeDependentBC_initialRateTimeHistory_scalar : pylith_fekernels_TimeDependentBC_initialRateTimeHistory_vector;
        break;
    case 0x0:
        PYLITH_JOURNAL_WARNING("Dirichlet BC provides no constraints.");
        break;
    default:
        PYLITH_JOURNAL_ERROR("Unknown combination of flags for Dirichlet BC terms (useInitial="<<_useInitial<<", useRate="<<_useRate<<", useTimeHistory="<<_useTimeHistory<<").");
        throw std::logic_error("Unknown combination of flags for Dirichlet BC terms.");
    } // switch

    void* context = NULL;
    const int labelId = 1;
    const int fieldIndex = solution.subfieldInfo(_field.c_str()).index;
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL_FIELD, _label.c_str(), _label.c_str(), fieldIndex, _constrainedDOF.size(), &_constrainedDOF[0], (void (*)())bcKernel, 1, &labelId, context); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsConstraint

// End of file
