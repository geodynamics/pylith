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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/materials/IsotropicLinearIncompElasticity.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactory
#include "pylith/fekernels/IsotropicLinearIncompElasticity.hh" // USES IsotropicLinearIncompElasticity kernels

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> \
    // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearIncompElasticity::IsotropicLinearIncompElasticity(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryElastic),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotopiclinearincomplinearelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearIncompElasticity::~IsotropicLinearIncompElasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearIncompElasticity::deallocate(void) {
    RheologyIncompressibleElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
void
pylith::materials::IsotropicLinearIncompElasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
bool
pylith::materials::IsotropicLinearIncompElasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicLinearIncompElasticity::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearIncompElasticity::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress();
        _auxiliaryFactory->addReferenceStrain();
    } // if
    _auxiliaryFactory->addShearModulus();
    _auxiliaryFactory->addBulkModulus();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearIncompElasticity::getKernelRHSResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointFunc g1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::g1u :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::g1u :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::g1u_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::g1u_refstate :
        NULL;

    PYLITH_METHOD_RETURN(g1u);
} // getKernelRHSResidualStress


// ---------------------------------------------------------------------------------------------------------------------
// Get pressure kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearIncompElasticity::getKernelRHSResidualPressure(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSResidualPressure(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointFunc g0p = (!_useReferenceState) ?
                         pylith::fekernels::IsotropicLinearIncompElasticity::g0p :
                         pylith::fekernels::IsotropicLinearIncompElasticity::g0p_refstate;

    PYLITH_METHOD_RETURN(g0p);
} // getKernelRHSREsidualPressure


// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearIncompElasticity::getKernelRHSJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointJac Jg3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::Jg3uu :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::Jg3uu :
        NULL;

    PYLITH_METHOD_RETURN(Jg3uu);
} // getKernelRHSJacobianElasticConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get inverse of the bulk modulus kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearIncompElasticity::getKernelRHSJacobianInverseBulkModulus(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSJacobianInverseBulkModulus(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg0pp = pylith::fekernels::IsotropicLinearIncompElasticity::Jg0pp;

    PYLITH_METHOD_RETURN(Jg0pp);
} // getKernelRHSJacobianInverseBulkModulus


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearIncompElasticity::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::stress :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::stress :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::stress_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::stress_refstate :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// End of file
