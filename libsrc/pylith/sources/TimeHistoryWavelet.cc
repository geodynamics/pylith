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
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "pylith/sources/AuxiliaryFactorySourceTime.hh" // USES AuxiliaryFactorySourceTime

#include "pylith/sources/TimeHistoryWavelet.hh" // implementation of object methods

#include "pylith/fekernels/TimeHistoryWavelet.hh" // USES TimeHistoryWavelet kernels
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::TimeHistoryWavelet::TimeHistoryWavelet(void) :
    _dbTimeHistory(NULL),
    _auxiliaryFactory(new pylith::sources::AuxiliaryFactorySourceTime), \
    _useTimeHistory(true) {
    pylith::utils::PyreComponent::setName("timehistorywavelet");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::TimeHistoryWavelet::~TimeHistoryWavelet(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::sources::TimeHistoryWavelet::deallocate(void) {
    SourceTimeFunctionMomentTensorForce::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set time history database.
void
pylith::sources::TimeHistoryWavelet::setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th) {
    PYLITH_COMPONENT_DEBUG("setTimeHistoryDB(th"<<th<<")");

    _dbTimeHistory = th;
} // setTimeHistoryDB


// ---------------------------------------------------------------------------------------------------------------------
// Get time history database.
const spatialdata::spatialdb::TimeHistory*
pylith::sources::TimeHistoryWavelet::getTimeHistoryDB(void) {
    return _dbTimeHistory;
} // getTimeHistoryDB


// ---------------------------------------------------------------------------------------------------------------------
// Use time history term in time history expression.
void
pylith::sources::TimeHistoryWavelet::useTimeHistory(const bool value) {
    PYLITH_COMPONENT_DEBUG("useTimeHistory(value="<<value<<")");

    _useTimeHistory = value;
} // useTimeHistory


// ---------------------------------------------------------------------------------------------------------------------
// Get flag associated with using time history term in time history expression.
bool
pylith::sources::TimeHistoryWavelet::useTimeHistory(void) const { // useTimeHistory
    return _useTimeHistory;
} // useTimeHistory


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::sources::AuxiliaryFactorySourceTime*
pylith::sources::TimeHistoryWavelet::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add source time subfields to auxiliary field.
void
pylith::sources::TimeHistoryWavelet::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    // _auxiliaryFactory->addTimeHistoryAmplitude(); // numA - 3
    _auxiliaryFactory->addTimeHistoryStartTime(); // numA - 2
    _auxiliaryFactory->addTimeHistoryValue(); // numA - 1
    if (_dbTimeHistory) {
        _dbTimeHistory->open();
    } // if
    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ---------------------------------------------------------------------------------------------------------------------
// Get g1v kernel for residual, G(t,s).
PetscPointFunc
pylith::sources::TimeHistoryWavelet::getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1v_explicit(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc g1v =
        (3 == spaceDim) ? pylith::fekernels::TimeHistoryWavelet3D::g1v :
        (2 == spaceDim) ? pylith::fekernels::TimeHistoryWaveletPlaneStrain::g1v :
        NULL;

    PYLITH_METHOD_RETURN(g1v);
} // getKernelResidualStress


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary fields.
void
pylith::sources::TimeHistoryWavelet::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                          const PylithReal t,
                                                          const PylithReal timeScale) {
    if (_useTimeHistory) {
        // assert(_normalizer);
        // const PylithScalar timeScale = _normalizer->getTimeScale();
        AuxiliaryFactorySourceTime::updateAuxiliaryField(auxiliaryField, t, timeScale, _dbTimeHistory);
    } // if
}


// End of file
