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

#include "NeumannTimeDependent.hh" // implementation of object methods

#include "TimeDependentAuxiliaryFactory.hh" // USES TimeDependentAuxiliaryFactory

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/fekernels/NeumannTimeDependent.hh" // USES NeumannTimeDependent kernels

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
const char* pylith::bc::NeumannTimeDependent::_pyreComponent = "neumanntimedependent";

// Local "private" functions.
namespace pylith {
    namespace bc {
        static void _setFEKernelsRHSResidualScalar(const NeumannTimeDependent* const bc,
                                                   PetscDS prob,
                                                   const PylithInt fieldIndex);

        static void _setFEKernelsRHSResidualVector(const NeumannTimeDependent* const bc,
                                                   PetscDS prob,
                                                   const PylithInt fieldIndex);
    } // bc
} // pylith

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::NeumannTimeDependent::NeumannTimeDependent(void) :
    _dbTimeHistory(NULL),
    _auxTimeDependentFactory(new pylith::bc::TimeDependentAuxiliaryFactory(pylith::bc::TimeDependentAuxiliaryFactory::TANGENTIAL_NORMAL)),
    _useInitial(true),
    _useRate(false),
    _useTimeHistory(false)
{ // constructor
    PyreComponent::name(_pyreComponent);
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::bc::NeumannTimeDependent::~NeumannTimeDependent(void)
{ // destructor
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::NeumannTimeDependent::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    Neumann::deallocate();
    delete _auxTimeDependentFactory; _auxTimeDependentFactory = NULL;

    _dbTimeHistory = NULL; // :KLUDGE: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set time history database.
void
pylith::bc::NeumannTimeDependent::dbTimeHistory(spatialdata::spatialdb::TimeHistory* th)
{ // dbTimeHistory
    _dbTimeHistory = th;
} // dbTimeHistory


// ----------------------------------------------------------------------
// Get time history database.
const spatialdata::spatialdb::TimeHistory*
pylith::bc::NeumannTimeDependent::dbTimeHistory(void)
{ // dbTimeHistory
    return _dbTimeHistory;
} // dbTimeHistory

// ----------------------------------------------------------------------
// Use initial value term in time history expression.
void
pylith::bc::NeumannTimeDependent::useInitial(const bool value)
{ // useInitial
    PYLITH_COMPONENT_DEBUG("useInitial(value="<<value<<")");

    _useInitial = value;
} // useInitial


// ----------------------------------------------------------------------
// Get flag associated with using initial value term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useInitial(void) const
{ // useInitial
    return _useInitial;
} // useInitial


// ----------------------------------------------------------------------
// Use rate value term in time history expression.
void
pylith::bc::NeumannTimeDependent::useRate(const bool value)
{ // useRate
    PYLITH_COMPONENT_DEBUG("useRate(value="<<value<<")");

    _useRate = value;
} // useRate


// ----------------------------------------------------------------------
// Get flag associated with using rate value term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useRate(void) const
{ // useRate
    return _useRate;
} // useRate


// ----------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::bc::NeumannTimeDependent::useTimeHistory(const bool value)
{ // useTimeHistory
    PYLITH_COMPONENT_DEBUG("useTimeHistory(value="<<value<<")");

    _useTimeHistory = value;
} // useTimeHistory


// ----------------------------------------------------------------------
// Get flag associated with using time history term in time history expression.
bool
pylith::bc::NeumannTimeDependent::useTimeHistory(void) const
{ // useTimeHistory
    return _useTimeHistory;
} // useTimeHistory


// ----------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::bc::NeumannTimeDependent::prestep(const double t,
                                          const double dt)
{ // prestep
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_prestep(t="<<t<<", dt="<<dt<<")");

    if (_useTimeHistory) {
        assert(_normalizer);
        assert(_auxField);

        const PylithScalar timeScale = _normalizer->timeScale();

        PetscErrorCode err = 0;

        PetscSection auxFieldsSection = _auxField->localSection(); assert(auxFieldsSection);
        PetscInt pStart = 0, pEnd = 0;
        err = PetscSectionGetChart(auxFieldsSection, &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
        pylith::topology::VecVisitorMesh auxFieldsVisitor(*_auxField);
        PetscScalar* auxFieldsArray = auxFieldsVisitor.localArray(); assert(auxFieldsArray);

        // Compute offset of time history subfields in auxiliary field.
        const PetscInt numComponents = _description.numComponents; assert(numComponents > 0);
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
pylith::bc::NeumannTimeDependent::_auxFieldSetup(const pylith::topology::Field& solution)
{ // _auxFieldSetup
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldSetup(solution="<<solution.label()<<")");

    assert(_auxTimeDependentFactory);
    assert(_normalizer);
    pylith::topology::Field::Description description = solution.subfieldInfo(_field.c_str()).description;
    if (_scaleName == std::string("pressure")) {
        description.scale = _normalizer->pressureScale();
    } else if (_scaleName == std::string("velocity")) {
        description.scale = _normalizer->lengthScale() / _normalizer->pressureScale();
    } else if (_scaleName == std::string("length")) {
        description.scale = _normalizer->pressureScale();
    } else if (_scaleName == std::string("time")) {
        description.scale = _normalizer->pressureScale();
    } else if (_scaleName == std::string("debsuty")) {
        description.scale = _normalizer->pressureScale();
    } else {
        std::ostringstream msg;
        msg << "Unknown name of scale ("<<_scaleName<<") for Neumann boundary condition '" << label() << "'.";
        PYLITH_COMPONENT_ERROR(msg.str());
        throw std::logic_error(msg.str());
    } // if/else
    _auxTimeDependentFactory->initialize(_auxField, *_normalizer, solution.spaceDim(), &description);

    // :ATTENTION: The order of the factory methods must match the order of the auxiliary subfields in the FE kernels.

    if (_useInitial) {
        _auxTimeDependentFactory->initialAmplitude();
    } // if
    if (_useRate) {
        _auxTimeDependentFactory->rateAmplitude();
        _auxTimeDependentFactory->rateStartTime();
    } // _useRate
    if (_useTimeHistory) {
        _auxTimeDependentFactory->timeHistoryAmplitude();
        _auxTimeDependentFactory->timeHistoryStartTime();
        _auxTimeDependentFactory->timeHistoryValue();
    } // _useTimeHistory

    PYLITH_METHOD_END;
}     // _auxFieldSetup


// ----------------------------------------------------------------------
// Get factory for setting up auxliary fields.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::NeumannTimeDependent::_auxFactory(void) {
    return _auxTimeDependentFactory;
} // _auxFactory


// ----------------------------------------------------------------------
// Does boundary conditon have point-wise functions (kernels) for integration/projection.
bool
pylith::bc::NeumannTimeDependent::_hasFEKernels(IntegratorPointwise::FEKernelKeys kernelsKey) const {
    bool hasKernels = false;
    switch (kernelsKey) {
    case KERNELS_RHS_RESIDUAL:
        hasKernels = true;
        break;
    case KERNELS_LHS_RESIDUAL:
    case KERNELS_RHS_JACOBIAN:
    case KERNELS_LHS_JACOBIAN:
    case KERNELS_LHS_JACOBIAN_LUMPEDINV:
    case KERNELS_UPDATE_STATE_VARS:
    case KERNELS_DERIVED_FIELDS:
        hasKernels = false;
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unrecognized finite-element kernels key '"<<kernelsKey<<"'.");
        throw std::logic_error("Unrecognized finite-element kernels key.");
    } // switch
    return hasKernels;
} // _hasKernels


// ----------------------------------------------------------------------
// Set point-wise functions (kernels) for integration/projection.
void
pylith::bc::NeumannTimeDependent::_setFEKernels(const pylith::topology::Field& solution,
                                                IntegratorPointwise::FEKernelKeys kernelsKey) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernels(solution="<<solution.label()<<", kernelsKey="<<kernelsKey<<")");

    const PetscDM dmSoln = solution.dmMesh(); assert(dmSoln);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    const bool isScalarField = _description.vectorFieldType == pylith::topology::Field::SCALAR;
    const int fieldIndex = solution.subfieldInfo(_field.c_str()).index;

    switch (kernelsKey) {
    case KERNELS_RHS_RESIDUAL:
        if (isScalarField) {
            _setFEKernelsRHSResidualScalar(this, prob, fieldIndex);
        } else {
            _setFEKernelsRHSResidualVector(this, prob, fieldIndex);
        } // if/else
        break;
    case KERNELS_LHS_RESIDUAL:
    case KERNELS_RHS_JACOBIAN:
    case KERNELS_LHS_JACOBIAN:
    case KERNELS_LHS_JACOBIAN_LUMPEDINV:
    case KERNELS_UPDATE_STATE_VARS:
    case KERNELS_DERIVED_FIELDS:
        break;
    } // switch
} // _setFEKernels

// ----------------------------------------------------------------------
// Set point-wise functions (kernels) for integrating RHS residual for scalar field.
void
pylith::bc::_setFEKernelsRHSResidualScalar(const NeumannTimeDependent* const bc,
                                           PetscDS prob,
                                           const PylithInt fieldIndex) {
    PYLITH_METHOD_BEGIN;

    PetscBdPointFunc g0 = NULL;
    PetscBdPointFunc g1 = NULL;

    const int bitInitial = bc->useInitial() ? 0x1 : 0x0;
    const int bitRate = bc->useRate() ? 0x2 : 0x0;
    const int bitTimeHistory = bc->useTimeHistory() ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;
    switch (bitUse) {
    case 0x1:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_initial_scalar;
        break;
    case 0x2:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_rate_scalar;
        break;
    case 0x4:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_timeHistory_scalar;
        break;
    case 0x3:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_initialRate_scalar;
        break;
    case 0x5:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_initialTimeHistory_scalar;
        break;
    case 0x6:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_rateTimeHistory_scalar;
        break;
    case 0x7:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_initialRateTimeHistory_scalar;
        break;
    case 0x0:
        //PYLITH_COMPONENT_WARNING("Neumann time-dependent BC provides no values.");
        break;
    default:
        throw std::logic_error("Unknown combination of flags for Neumann time-dependent BC terms.");
    } // switch

    PetscErrorCode err = PetscDSSetBdResidual(prob, fieldIndex, g0, g1); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidualScalar


// ----------------------------------------------------------------------
// Set point-wise functions (kernels) for integrating RHS residual for vector field.
void
pylith::bc::_setFEKernelsRHSResidualVector(const NeumannTimeDependent* const bc,
                                           PetscDS prob,
                                           const PylithInt fieldIndex) {
    PYLITH_METHOD_BEGIN;

    PetscBdPointFunc g0 = NULL;
    PetscBdPointFunc g1 = NULL;

    const int bitInitial = bc->useInitial() ? 0x1 : 0x0;
    const int bitRate = bc->useRate() ? 0x2 : 0x0;
    const int bitTimeHistory = bc->useTimeHistory() ? 0x4 : 0x0;
    const int bitUse = bitInitial | bitRate | bitTimeHistory;
    switch (bitUse) {
    case 0x1:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_initial_vector;
        break;
    case 0x2:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_rate_vector;
        break;
    case 0x4:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_timeHistory_vector;
        break;
    case 0x3:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_initialRate_vector;
        break;
    case 0x5:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_initialTimeHistory_vector;
        break;
    case 0x6:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_rateTimeHistory_vector;
        break;
    case 0x7:
        g0 = pylith::fekernels::NeumannTimeDependent::g0_initialRateTimeHistory_vector;
        break;
    case 0x0:
        //PYLITH_COMPONENT_WARNING("Neumann time-dependent BC provides no values.");
        break;
    default:
        throw std::logic_error("Unknown combination of flags for Neumann time-dependent BC terms.");
    } // switch

    PetscErrorCode err = PetscDSSetBdResidual(prob, fieldIndex, g0, g1); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidualVector


// End of file
