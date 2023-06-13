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

#include "pylith/sources/GaussianWavelet.hh" // implementation of object methods

#include "pylith/sources/AuxiliaryFactorySourceTime.hh" // USES AuxiliaryFactorySourceTime
#include "pylith/fekernels/GaussianWavelet.hh" // USES GaussianWavelet kernels
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::GaussianWavelet::GaussianWavelet(void) :
    _auxiliaryFactory(new pylith::sources::AuxiliaryFactorySourceTime) {
    pylith::utils::PyreComponent::setName("gaussianwavelet");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::GaussianWavelet::~GaussianWavelet(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::sources::GaussianWavelet::deallocate(void) {
    SourceTimeFunctionMomentTensorForce::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::sources::AuxiliaryFactoryMomentTensorForce*
pylith::sources::GaussianWavelet::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add source time subfields to auxiliary field.
void
pylith::sources::GaussianWavelet::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    _auxiliaryFactory->addCenterFrequency(); // numA - 1

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ---------------------------------------------------------------------------------------------------------------------
// Get g1v kernel for residual, G(t,s).
PetscPointFunc
pylith::sources::GaussianWavelet::getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1v_explicit(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc g1v =
        (3 == spaceDim) ? pylith::fekernels::GaussianWavelet3D::g1v :
        (2 == spaceDim) ? pylith::fekernels::GaussianWaveletPlaneStrain::g1v :
        NULL;

    PYLITH_METHOD_RETURN(g1v);
} // getKernelResidualStress


// End of file
