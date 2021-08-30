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

#include "pylith/materials/IsotropicPowerLaw.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryViscoelastic.hh" // USES AuxiliaryFactoryViscoelastic
#include "pylith/fekernels/IsotropicPowerLaw.hh" // USES IsotropicPowerLaw kernels
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicPowerLaw::IsotropicPowerLaw(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryViscoelastic),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotropicpowerlaw");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicPowerLaw::~IsotropicPowerLaw(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicPowerLaw::deallocate(void) {
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicPowerLaw::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicPowerLaw::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicPowerLaw::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicPowerLaw::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    _auxiliaryFactory->addShearModulus();
    _auxiliaryFactory->addBulkModulus();
    _auxiliaryFactory->addPowerLawReferenceStrainRate();
    _auxiliaryFactory->addPowerLawReferenceStress();
    _auxiliaryFactory->addPowerLawExponent();
    _auxiliaryFactory->addViscousStrain();
    _auxiliaryFactory->addStress();
    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress();
        _auxiliaryFactory->addReferenceStrain();
    } // if

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFunc
pylith::materials::IsotropicPowerLaw::getKernelResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::f1v :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::f1v :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::f1v_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::f1v_refstate :
        NULL;

    PYLITH_METHOD_RETURN(f1u);
} // getKernelResidualStress


// ------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJac
pylith::materials::IsotropicPowerLaw::getKernelJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3uu =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::Jf3vu :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::Jf3vu :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::Jf3vu_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::Jf3vu_refstate :
        NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelJacobianElasticConstants


// ------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicPowerLaw::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::cauchyStress :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::cauchyStress :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::cauchyStress_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::cauchyStress_refstate :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::IsotropicPowerLaw::updateKernelConstants(pylith::real_array* kernelConstants,
                                                            const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants"<<kernelConstants<<", dt="<<dt<<")");
    // ******** Should alpha (time integration parameter) be included here?

    assert(kernelConstants);

    if (1 != kernelConstants->size()) { kernelConstants->resize(1);}
    (*kernelConstants)[0] = dt;

    PYLITH_METHOD_END;
} // updateKernelConstants


// ------------------------------------------------------------------------------------------------
// Add kernels for updating state variables.
void
pylith::materials::IsotropicPowerLaw::addKernelsUpdateStateVars(std::vector<ProjectKernels>* kernels,
                                                                const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVars(kernels="<<kernels<<", coordsys="<<coordsys<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc funcViscousStrain =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::updateViscousStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::updateViscousStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::updateViscousStrain_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::updateViscousStrain_refstate :
        NULL;
    const PetscPointFunc funcStress =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::updateStress :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::updateStress :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicPowerLaw3D::updateStress_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicPowerLawPlaneStrain::updateStress_refstate :
        NULL;

    assert(kernels);
    size_t prevNumKernels = kernels->size();
    kernels->resize(prevNumKernels + 2);
    (*kernels)[prevNumKernels+0] = ProjectKernels("viscous_strain", funcViscousStrain);
    (*kernels)[prevNumKernels+1] = ProjectKernels("stress", funcStress);

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVars


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFunc
pylith::materials::IsotropicPowerLaw::getInterfaceKernelResidualF0Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelResidualF0Neg().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFunc
pylith::materials::IsotropicPowerLaw::getInterfaceKernelResidualF0Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelResidualF0Pos().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get f1 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFunc
pylith::materials::IsotropicPowerLaw::getInterfaceKernelResidualF1Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernel10Neg().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get f1 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFunc
pylith::materials::IsotropicPowerLaw::getInterfaceKernelResidualF1Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelResidualF1Pos().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get Jf1lu kernel for LHS Jacobian F(t,s,dot{s}) for negative fault face.
PetscBdPointJac
pylith::materials::IsotropicPowerLaw::getInterfaceKernelJacobianF1Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelJacobianF1Neg().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get Jf1lu kernel for LHS Jacobian F(t,s,dot{s}) for positive fault face.
PetscBdPointJac
pylith::materials::IsotropicPowerLaw::getInterfaceKernelJacobianF1Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelJacobianF1Pos().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get Jf3lu kernel for LHS Jacobian F(t,s,dot{s}) for negative fault face.
PetscBdPointJac
pylith::materials::IsotropicPowerLaw::getInterfaceKernelJacobianF3Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelJacobianF3Neg().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get Jf3lu kernel for LHS Jacobian F(t,s,dot{s}) for positive fault face.
PetscBdPointJac
pylith::materials::IsotropicPowerLaw::getInterfaceKernelJacobianF3Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelJacobianF3Pos().");
    return NULL;
}


// End of file
