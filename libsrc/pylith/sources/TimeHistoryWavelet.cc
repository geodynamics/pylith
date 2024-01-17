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

#include "pylith/sources/TimeHistoryWavelet.hh" // implementation of object methods

#include "pylith/sources/AuxiliaryFactorySourceTime.hh" // USES AuxiliaryFactorySourceTime
#include "pylith/fekernels/TimeHistoryWavelet.hh" // USES TimeHistoryWavelet kernels
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::TimeHistoryWavelet::TimeHistoryWavelet(void) :
    _auxiliaryFactory(new pylith::sources::AuxiliaryFactorySourceTime) {
    pylith::utils::PyreComponent::setName("timehistorywavelet");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::TimeHistoryWavelet::~TimeHistoryWavelet(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Set time history database.
void
pylith::sources::SquarePulseSource::setTimeHistoryDB(spatialdata::spatialdb::TimeHistory *th) {
    PYLITH_COMPONENT_DEBUG("setTimeHistoryDB(th" << th << ")");

    _dbTimeHistory = th;
} // setTimeHistoryDB


// ---------------------------------------------------------------------------------------------------------------------
// Get time history database.
const spatialdata::spatialdb::TimeHistory *
pylith::sources::SquarePulseSource::getTimeHistoryDB(void) {
    return _dbTimeHistory;
} // getTimeHistoryDB


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::sources::TimeHistoryWavelet::deallocate(void) {
    SourceTimeFunctionMomentTensorForce::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::sources::AuxiliaryFactoryMomentTensorForce*
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

    _auxiliaryFactory->addCenterFrequency(); // numA - 1
    _auxiliaryFactory->addTimeHistoryAmplitude();
    _auxiliaryFactory->addTimeHistoryStartTime();
    _auxiliaryFactory->addTimeHistoryValue();
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


// End of file
