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

#include "pylith/materials/IsotropicLinearElasticity.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactoryElastic
#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity kernels
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearElasticity::IsotropicLinearElasticity(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryElastic),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotropiclinearelasticity");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearElasticity::~IsotropicLinearElasticity(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearElasticity::deallocate(void) {
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearElasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearElasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicLinearElasticity::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
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


// ------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFunc
pylith::materials::IsotropicLinearElasticity::getKernelResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1v :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1v :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1v_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1v_refstate :
        NULL;

    PYLITH_METHOD_RETURN(f1u);
} // getKernelResidualStress


// ------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearElasticity::getKernelJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::Jf3vu :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf3vu :
        NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelJacobianElasticConstants


// ------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearElasticity::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress_refstate :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearElasticity::getInterfaceKernelResidualF0Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelResidualF0Neg(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_neg :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_neg :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_refstate_neg :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_refstate_neg :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearElasticity::getInterfaceKernelResidualF0Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelResidualF0Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_pos :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_pos :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_refstate_pos :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_refstate_pos :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get f1 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearElasticity::getInterfaceKernelResidualF1Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelResidualF1Neg(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1l_neg :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1l_neg :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1l_refstate_neg :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1l_refstate_neg :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get f1 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearElasticity::getInterfaceKernelResidualF1Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelResidualF1Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1l_pos :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1l_pos :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1l_refstate_pos :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1l_refstate_pos :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get Jf1lu kernel for LHS Jacobian F(t,s,dot{s}) for negative fault face.
PetscBdPointJac
pylith::materials::IsotropicLinearElasticity::getInterfaceKernelJacobianF1Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelResidualF1Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointJac kernel =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::Jf1lu_neg :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf1lu_neg :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get Jf1lu kernel for LHS Jacobian F(t,s,dot{s}) for positive fault face.
PetscBdPointJac
pylith::materials::IsotropicLinearElasticity::getInterfaceKernelJacobianF1Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelResidualF1Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointJac kernel =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::Jf1lu_pos :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf1lu_pos :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get Jf3lu kernel for LHS Jacobian F(t,s,dot{s}) for negative fault face.
PetscBdPointJac
pylith::materials::IsotropicLinearElasticity::getInterfaceKernelJacobianF3Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelResidualF1Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointJac kernel =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::Jf3lu_neg :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf3lu_neg :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get Jf3lu kernel for LHS Jacobian F(t,s,dot{s}) for positive fault face.
PetscBdPointJac
pylith::materials::IsotropicLinearElasticity::getInterfaceKernelJacobianF3Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelResidualF1Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointJac kernel =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::Jf3lu_pos :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf3lu_pos :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// End of file
