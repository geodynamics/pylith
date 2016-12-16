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
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
const char* pylith::bc::DirichletTimeDependent::_pyreComponent = "dirichlettimedependent";

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::DirichletTimeDependent::DirichletTimeDependent(void) :
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

    PYLITH_METHOD_END;
} // deallocate


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
    const size_t numDOFConstrained = _constrainedDOF.size();

    // :ATTENTION: The order for subfieldAdd() must match the order of the auxiliary fields in the FE kernels.

    // Initial amplitude
    if (_useInitial) {
        const char* allCompNames[3] = {"initial_amplitude_x", "initial_amplitude_y", "initial_amplitude_z"};
        const char** constrainedCompNames = (numDOFConstrained) ? new const char*[numDOFConstrained] : NULL;
        for (size_t i=0; i < numDOFConstrained; ++i) {
            constrainedCompNames[i] = allCompNames[_constrainedDOF[i]];
        } // for
        const pylith::topology::Field::DiscretizeInfo& initialAmplitudeFEInfo = this->auxFieldDiscretization("initial_amplitude");
        _auxFields->subfieldAdd("initial_amplitude", constrainedCompNames, numDOFConstrained, pylith::topology::Field::OTHER, initialAmplitudeFEInfo.basisOrder, initialAmplitudeFEInfo.quadOrder, initialAmplitudeFEInfo.isBasisContinuous, lengthScale);
        _auxFieldsQuery->queryFn("initial_amplitude", pylith::topology::FieldQuery::dbQueryGeneric);
        delete[] constrainedCompNames; constrainedCompNames = NULL;
    } // if


    // Rate amplitude and start time.
    if (_useRate) {
        const char* allCompNames[3] = {"rate_amplitude_x", "rate_amplitude_y", "rate_amplitude_z"};
        const char** constrainedCompNames = (numDOFConstrained) ? new const char*[numDOFConstrained] : NULL;
        for (size_t i=0; i < numDOFConstrained; ++i) {
            constrainedCompNames[i] = allCompNames[_constrainedDOF[i]];
        } // for
        const pylith::topology::Field::DiscretizeInfo& rateAmplitudeFEInfo = this->auxFieldDiscretization("rate_amplitude");
        _auxFields->subfieldAdd("rate_amplitude", constrainedCompNames, numDOFConstrained, pylith::topology::Field::OTHER, rateAmplitudeFEInfo.basisOrder, rateAmplitudeFEInfo.quadOrder, rateAmplitudeFEInfo.isBasisContinuous, velocityScale);
        _auxFieldsQuery->queryFn("rate_amplitude", pylith::topology::FieldQuery::dbQueryGeneric);
        delete[] constrainedCompNames; constrainedCompNames = NULL;

        const char* rateStartNames[1] = {"rate_start"};
        const pylith::topology::Field::DiscretizeInfo& rateStartFEInfo = this->auxFieldDiscretization("rate_start");
        _auxFields->subfieldAdd("rate_start", rateStartNames, 1, pylith::topology::Field::SCALAR, rateStartFEInfo.basisOrder, rateStartFEInfo.quadOrder, rateStartFEInfo.isBasisContinuous, timeScale);
        _auxFieldsQuery->queryFn("rate_start", pylith::topology::FieldQuery::dbQueryGeneric);
    } // if


    // Time history amplitude and start time.
    if (_useTimeHistory) {
        const char* allCompNames[3] = {"time_history_amplitude_x", "time_history_amplitude_y", "time_history_amplitude_z"};
        const char** constrainedCompNames = (numDOFConstrained) ? new const char*[numDOFConstrained] : NULL;
        for (size_t i=0; i < numDOFConstrained; ++i) {
            constrainedCompNames[i] = allCompNames[_constrainedDOF[i]];
        } // for
        const pylith::topology::Field::DiscretizeInfo& timeHistoryAmplitudeFEInfo = this->auxFieldDiscretization("time_history_amplitude");
        _auxFields->subfieldAdd("time_history_amplitude", constrainedCompNames, numDOFConstrained, pylith::topology::Field::OTHER, timeHistoryAmplitudeFEInfo.basisOrder, timeHistoryAmplitudeFEInfo.quadOrder, timeHistoryAmplitudeFEInfo.isBasisContinuous, lengthScale);
        _auxFieldsQuery->queryFn("time_history_amplitude", pylith::topology::FieldQuery::dbQueryGeneric);
        delete[] constrainedCompNames; constrainedCompNames = NULL;

        const char* timeHistoryStartNames[1] = {"time_history_start"};
        const pylith::topology::Field::DiscretizeInfo& timeHistoryStartFEInfo = this->auxFieldDiscretization("time_history_start");
        _auxFields->subfieldAdd("time_history_start", timeHistoryStartNames, 1, pylith::topology::Field::SCALAR, timeHistoryStartFEInfo.basisOrder, timeHistoryStartFEInfo.quadOrder, timeHistoryStartFEInfo.isBasisContinuous, timeScale);
        _auxFieldsQuery->queryFn("time_history_start", pylith::topology::FieldQuery::dbQueryGeneric);

        const char* timeHistoryValueNames[1] = {"time_history_value"};
        _auxFields->subfieldAdd("time_history_value", timeHistoryValueNames, 1, pylith::topology::Field::SCALAR, timeHistoryAmplitudeFEInfo.basisOrder, timeHistoryAmplitudeFEInfo.quadOrder, timeHistoryAmplitudeFEInfo.isBasisContinuous, 1.0);
        _auxFieldsQuery->queryFn("time_history_value", NULL);
    } // if


    PYLITH_METHOD_END;
}     // _auxFieldsSetup


// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::bc::DirichletTimeDependent::_setFEKernelsConstraint(const topology::Field& solution)
{ // _setFEKernelsConstraint
    PYLITH_JOURNAL_DEBUG("setFEKernelsConstraint(solution="<<solution.label()<<")");

    PYLITH_JOURNAL_ERROR(":TODO: @brad Implement _setFEKernelsConstraint().");
} // _setFEKernelsConstraint

// End of file
