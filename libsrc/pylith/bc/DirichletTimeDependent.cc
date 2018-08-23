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

#include "TimeDependentAuxiliaryFactory.hh" // USES TimeDependentAuxiliaryFactory

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/fekernels/TimeDependentFn.hh" // USES TimeDependentFn kernels

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
const char* pylith::bc::DirichletTimeDependent::_pyreComponent = "dirichlettimedependent";

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletTimeDependent::DirichletTimeDependent(void) :
    _dbTimeHistory(NULL),
    _auxTimeDependentFactory(new pylith::bc::TimeDependentAuxiliaryFactory),
    _bcKernel(NULL),
    _useInitial(true),
    _useRate(false),
    _useTimeHistory(false) { // constructor
    PyreComponent::name(_pyreComponent);
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::bc::DirichletTimeDependent::~DirichletTimeDependent(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::DirichletTimeDependent::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    ConstraintBoundary::deallocate();

    _dbTimeHistory = NULL; // :KLUDGE: Use shared pointer.
    delete _auxTimeDependentFactory;_auxTimeDependentFactory = NULL;
    _bcKernel = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set time history database.
void
pylith::bc::DirichletTimeDependent::dbTimeHistory(spatialdata::spatialdb::TimeHistory* th) {
    _dbTimeHistory = th;
} // dbTimeHistory


// ----------------------------------------------------------------------
// Get time history database.
const spatialdata::spatialdb::TimeHistory*
pylith::bc::DirichletTimeDependent::dbTimeHistory(void) {
    return _dbTimeHistory;
} // dbTimeHistory


// ----------------------------------------------------------------------
// Use initial value term in time history expression.
void
pylith::bc::DirichletTimeDependent::useInitial(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInitial(value="<<value<<")");

    _useInitial = value;
} // useInitial


// ----------------------------------------------------------------------
// Get flag associated with using initial value term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useInitial(void) const {
    return _useInitial;
} // useInitial


// ----------------------------------------------------------------------
// Use rate value term in time history expression.
void
pylith::bc::DirichletTimeDependent::useRate(const bool value) {
    PYLITH_COMPONENT_DEBUG("useRate(value="<<value<<")");

    _useRate = value;
} // useRate


// ----------------------------------------------------------------------
// Get flag associated with using rate value term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useRate(void) const {
    return _useRate;
} // useRate


// ----------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::bc::DirichletTimeDependent::useTimeHistory(const bool value) {
    PYLITH_COMPONENT_DEBUG("useTimeHistory(value="<<value<<")");

    _useTimeHistory = value;
} // useTimeHistory


// ----------------------------------------------------------------------
// Get flag associated with using time history term in time history expression.
bool
pylith::bc::DirichletTimeDependent::useTimeHistory(void) const {
    return _useTimeHistory;
} // useTimeHistory


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::bc::DirichletTimeDependent::prestep(const double t,
                                            const double dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_prestep(t="<<t<<", dt="<<dt<<")");

    if (_useTimeHistory) {
        assert(_normalizer);
        assert(_auxField);

        const PylithScalar timeScale = _normalizer->timeScale();

        PetscErrorCode err;

        PetscSection auxFieldsSection = _auxField->localSection();assert(auxFieldsSection);
        PetscInt pStart = 0, pEnd = 0;
        err = PetscSectionGetChart(auxFieldsSection, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
        pylith::topology::VecVisitorMesh auxFieldsVisitor(*_auxField);
        PetscScalar* auxFieldsArray = auxFieldsVisitor.localArray();assert(auxFieldsArray);

        const PylithInt numComponents = _description.numComponents;assert(numComponents > 0);
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

    //_auxField->view("AUXILIARY FIELD"); // :DEBUG: TEMPORARY

    PYLITH_METHOD_END;
} // prestep


// ----------------------------------------------------------------------
// Get factory for setting up auxliary fields.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::DirichletTimeDependent::_auxFactory(void) {
    return _auxTimeDependentFactory;
} // auxFactory


// ----------------------------------------------------------------------
// Setup auxiliary subfields (discretization and query fns).
void
pylith::bc::DirichletTimeDependent::_auxFieldSetup(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldsSetup(solution="<<solution.label()<<")");

    assert(_auxTimeDependentFactory);
    assert(_normalizer);
    _auxTimeDependentFactory->initialize(_auxField, *_normalizer, solution.spaceDim(),
                                         &solution.subfieldInfo(_field.c_str()).description);

    // :ATTENTION: The order of the factory methods must match the order of the auxiliary subfields in the FE kernels.

    if (_useInitial) {
        _auxTimeDependentFactory->addInitialAmplitude();
    } // if
    if (_useRate) {
        _auxTimeDependentFactory->addRateAmplitude();
        _auxTimeDependentFactory->addRateStartTime();
    } // _useRate
    if (_useTimeHistory) {
        _auxTimeDependentFactory->addTimeHistoryAmplitude();
        _auxTimeDependentFactory->addTimeHistoryStartTime();
        _auxTimeDependentFactory->addTimeHistoryValue();
    } // _useTimeHistory

    PYLITH_METHOD_END;
} // _auxFieldSetup


// ----------------------------------------------------------------------
// Set kernels for setting Dirhclet values.
void
pylith::bc::DirichletTimeDependent::_setFEKernelConstraint(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelConstraint(solution="<<solution.label()<<")");

    const PetscDM dmSoln = solution.dmMesh();assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);

    const bool isScalarField = _description.vectorFieldType == pylith::topology::Field::SCALAR;

    const int bitInitial = _useInitial ? 0x1 : 0x0;
    const int bitRate = _useRate ? 0x2 : 0x0;
    const int bitTimeHistory = _useTimeHistory ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;
    switch (bitUse) {
    case 0x1:
        _bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::initial_scalar : pylith::fekernels::TimeDependentFn::initial_vector;
        break;
    case 0x2:
        _bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::rate_scalar : pylith::fekernels::TimeDependentFn::rate_vector;
        break;
    case 0x4:
        _bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::timeHistory_scalar : pylith::fekernels::TimeDependentFn::timeHistory_vector;
        break;
    case 0x3:
        _bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::initialRate_scalar : pylith::fekernels::TimeDependentFn::initialRate_vector;
        break;
    case 0x5:
        _bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::initialTimeHistory_scalar : pylith::fekernels::TimeDependentFn::initialTimeHistory_vector;
        break;
    case 0x6:
        _bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::rateTimeHistory_scalar : pylith::fekernels::TimeDependentFn::rateTimeHistory_vector;
        break;
    case 0x7:
        _bcKernel = (isScalarField) ? pylith::fekernels::TimeDependentFn::initialRateTimeHistory_scalar : pylith::fekernels::TimeDependentFn::initialRateTimeHistory_vector;
        break;
    case 0x0:
        PYLITH_COMPONENT_WARNING("Dirichlet BC provides no constraints.");
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown combination of flags for Dirichlet BC terms (useInitial="<<_useInitial<<", useRate="<<_useRate<<", useTimeHistory="<<_useTimeHistory<<").");
        throw std::logic_error("Unknown combination of flags for Dirichlet BC terms.");
    } // switch

    PYLITH_METHOD_END;
} // _setFEKernelConstraint


// ----------------------------------------------------------------------
// Get point-wise function (kernel) for settings constraint from auxiliary field.
PetscPointFunc
pylith::bc::DirichletTimeDependent::_getFEKernelConstraint(void) {
    return _bcKernel;
} // _getFEKernelConstraint


// End of file
