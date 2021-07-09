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
pylith::materials::IsotropicLinearIncompElasticity::getKernelResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::f1u :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::f1u :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::f1u_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::f1u_refstate :
        NULL;

    PYLITH_METHOD_RETURN(f1u);
} // getKernelRHSResidualStress


// ---------------------------------------------------------------------------------------------------------------------
// Get pressure kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearIncompElasticity::getKernelResidualPressure(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSResidualPressure(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointFunc f0p = (!_useReferenceState) ?
                         pylith::fekernels::IsotropicLinearIncompElasticity::f0p :
                         pylith::fekernels::IsotropicLinearIncompElasticity::f0p_refstate;

    PYLITH_METHOD_RETURN(f0p);
} // getKernelResidualPressure


// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearIncompElasticity::getKernelJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jg3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::Jf3uu :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::Jf3uu :
        NULL;

    PYLITH_METHOD_RETURN(Jg3uu);
} // getKernelJacobianElasticConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get inverse of the bulk modulus kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearIncompElasticity::getKernelJacobianInverseBulkModulus(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJacobianInverseBulkModulus(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jf0pp = pylith::fekernels::IsotropicLinearIncompElasticity::Jf0pp;

    PYLITH_METHOD_RETURN(Jf0pp);
} // getKernelJacobianInverseBulkModulus


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearIncompElasticity::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::cauchyStress :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::cauchyStress :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticity3D::cauchyStress_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain::cauchyStress_refstate :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// End of file
