// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "KinSrcTimeHistory.hh" // implementation of object methods

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // USES KinSrcAuxiliaryFactory

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory#include <cassert> // USES assert()
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcTimeHistory::KinSrcTimeHistory(void) :
    _dbTimeHistory(NULL) {
    pylith::utils::PyreComponent::setName("kinsrctimehistory");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcTimeHistory::~KinSrcTimeHistory(void) {
    _dbTimeHistory = NULL; // :KLUDGE: Use shared pointer.
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set time history database.
void
pylith::faults::KinSrcTimeHistory::setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th) {
    PYLITH_COMPONENT_DEBUG("setTimeHistoryDB(th"<<th<<")");

    _dbTimeHistory = th;
} // setTimeHistoryDB


// ---------------------------------------------------------------------------------------------------------------------
// Get time history database.
const spatialdata::spatialdb::TimeHistory*
pylith::faults::KinSrcTimeHistory::getTimeHistoryDB(void) {
    return _dbTimeHistory;
} // getTimeHistoryDB


// ----------------------------------------------------------------------
// Get requested slip subfields at time t.
void
pylith::faults::KinSrcTimeHistory::getSlipSubfields(PetscVec slipLocalVec,
                                                    pylith::topology::Field* faultAuxiliaryField,
                                                    const PylithReal t,
                                                    const PylithReal dt,
                                                    const PylithReal timeScale,
                                                    const int bitSlipSubfields) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getSlipSubfields="<<slipLocalVec<<", faultAuxiliaryField="<<faultAuxiliaryField
                                              <<", t="<<t<<", dt="<<dt<<", timeScale="<<timeScale
                                              <<", bitSlipSubfields="<<bitSlipSubfields<<")");
    KinSrcAuxiliaryFactory::updateTimeHistoryValue(_auxiliaryField, t, timeScale, _dbTimeHistory);
    KinSrc::getSlipSubfields(slipLocalVec, faultAuxiliaryField, t, dt, timeScale, bitSlipSubfields);

    PYLITH_METHOD_END;
} // getSlipSubfields


// ---------------------------------------------------------------------------------------------------------------------
// Slip time function kernel.
void
pylith::faults::KinSrcTimeHistory::slipFn(const PylithInt dim,
                                          const PylithInt numS,
                                          const PylithInt numA,
                                          const PylithInt sOff[],
                                          const PylithInt sOff_x[],
                                          const PylithScalar s[],
                                          const PylithScalar s_t[],
                                          const PylithScalar s_x[],
                                          const PylithInt aOff[],
                                          const PylithInt aOff_x[],
                                          const PylithScalar a[],
                                          const PylithScalar a_t[],
                                          const PylithScalar a_x[],
                                          const PylithReal t,
                                          const PylithScalar x[],
                                          const PylithInt numConstants,
                                          const PylithScalar constants[],
                                          PylithScalar slip[]) {
    const PylithInt _numA = 3;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slip);

    const PylithInt i_initiationTime = 0;
    const PylithInt i_finalSlip = 1;
    const PylithInt i_timeHistoryValue = 2;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar timeHistoryValue = a[aOff[i_timeHistoryValue]];

    const PylithInt i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        for (PylithInt i = 0; i < dim; ++i) {
            slip[i] = finalSlip[i] * timeHistoryValue;
        } // for
    } // if

} // slipFn


// ---------------------------------------------------------------------------------------------------------------------
// Slip rate time function kernel.
void
pylith::faults::KinSrcTimeHistory::slipRateFn(const PylithInt dim,
                                              const PylithInt numS,
                                              const PylithInt numA,
                                              const PylithInt sOff[],
                                              const PylithInt sOff_x[],
                                              const PylithScalar s[],
                                              const PylithScalar s_t[],
                                              const PylithScalar s_x[],
                                              const PylithInt aOff[],
                                              const PylithInt aOff_x[],
                                              const PylithScalar a[],
                                              const PylithScalar a_t[],
                                              const PylithScalar a_x[],
                                              const PylithReal t,
                                              const PylithScalar x[],
                                              const PylithInt numConstants,
                                              const PylithScalar constants[],
                                              PylithScalar slipRate[]) {
    const PylithInt _numA = 3;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slipRate);

    const PylithInt i_initiationTime = 0;
    const PylithInt i_finalSlip = 1;
    const PylithInt i_timeHistoryValue = 2;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar timeHistoryValue = a[aOff[i_timeHistoryValue]];

    const PylithInt i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        for (PylithInt i = 0; i < dim; ++i) {
            slipRate[i] = finalSlip[i] * timeHistoryValue;
        } // for
    } // if

} // slipRateFn


// ---------------------------------------------------------------------------------------------------------------------
// Slip acceleration time function kernel.
void
pylith::faults::KinSrcTimeHistory::slipAccFn(const PylithInt dim,
                                             const PylithInt numS,
                                             const PylithInt numA,
                                             const PylithInt sOff[],
                                             const PylithInt sOff_x[],
                                             const PylithScalar s[],
                                             const PylithScalar s_t[],
                                             const PylithScalar s_x[],
                                             const PylithInt aOff[],
                                             const PylithInt aOff_x[],
                                             const PylithScalar a[],
                                             const PylithScalar a_t[],
                                             const PylithScalar a_x[],
                                             const PylithReal t,
                                             const PylithScalar x[],
                                             const PylithInt numConstants,
                                             const PylithScalar constants[],
                                             PylithScalar slipAcc[]) {
    const PylithInt _numA = 3;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slipAcc);

    const PylithInt i_initiationTime = 0;
    const PylithInt i_finalSlip = 1;
    const PylithInt i_timeHistoryValue = 2;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar timeHistoryValue = a[aOff[i_timeHistoryValue]];

    const PylithInt i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        for (PylithInt i = 0; i < dim; ++i) {
            slipAcc[i] = finalSlip[i] * timeHistoryValue;
        } // for
    } // if

} // slipAccFn


// ---------------------------------------------------------------------------------------------------------------------
// Preinitialize earthquake source. Set names/sizes of auxiliary subfields.
void
pylith::faults::KinSrcTimeHistory::_auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                                        const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxiliaryFieldSetup()");

    assert(_auxiliaryFactory);
    assert(cs);
    _auxiliaryFactory->initialize(_auxiliaryField, normalizer, cs->getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the slip time function
    // kernel.

    _auxiliaryFactory->addInitiationTime(); // 0
    _auxiliaryFactory->addFinalSlip(); // 1
    _auxiliaryFactory->addTimeHistoryValue(); // 2

    assert(_dbTimeHistory);
    _dbTimeHistory->open();

    _slipFnKernel = pylith::faults::KinSrcTimeHistory::slipFn;
    _slipRateFnKernel = pylith::faults::KinSrcTimeHistory::slipRateFn;
    _slipAccFnKernel = pylith::faults::KinSrcTimeHistory::slipAccFn;

    PYLITH_METHOD_END;
} // _auxiliaryFieldSetup


// End of file
