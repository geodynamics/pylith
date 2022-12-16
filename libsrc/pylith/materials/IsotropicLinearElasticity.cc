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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/materials/IsotropicLinearElasticity.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactoryElastic
#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity kernels
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearElasticity::IsotropicLinearElasticity(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryElastic),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotropiclinearelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearElasticity::~IsotropicLinearElasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearElasticity::deallocate(void) {
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearElasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearElasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicLinearElasticity::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearElasticity::addAuxiliarySubfields(void) {
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
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFunc
pylith::materials::IsotropicLinearElasticity::getKernelf1v(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f1v =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1v_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1v_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1v_infinitesimalStrain_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1v_infinitesimalStrain_refstate :
        NULL;

    PYLITH_METHOD_RETURN(f1v);
} // getKernelf1v


// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearElasticity::getKernelJf3vu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3vu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::Jf3vu_infinitesimalStrain :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf3vu_infinitesimalStrain :
        NULL;

    PYLITH_METHOD_RETURN(Jf3vu);
} // getKernelJacobianElasticConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearElasticity::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress_infinitesimalStrain_asVector :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress_infinitesimalStrain_asVector :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress_infinitesimalStrain_refstate_asVector :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress_infinitesimalStrain_refstate_asVector :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// End of file
