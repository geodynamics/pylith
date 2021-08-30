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

#include "pylith/materials/IsotropicLinearGenMaxwell.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryViscoelastic.hh" // USES AuxiliaryFactoryViscoelastic
#include "pylith/fekernels/IsotropicLinearGenMaxwell.hh" // USES IsotropicLinearGenMaxwell kernels
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearGenMaxwell::IsotropicLinearGenMaxwell(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryViscoelastic),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotropiclineargenmaxwell");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearGenMaxwell::~IsotropicLinearGenMaxwell(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearGenMaxwell::deallocate(void) {
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearGenMaxwell::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearGenMaxwell::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicLinearGenMaxwell::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearGenMaxwell::addAuxiliarySubfields(void) {
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
    _auxiliaryFactory->addMaxwellTimeGeneralizedMaxwell(); // 3
    _auxiliaryFactory->addShearModulusRatioGeneralizedMaxwell(); // 4
    _auxiliaryFactory->addViscousStrainGeneralizedMaxwell(); // 5
    _auxiliaryFactory->addTotalStrain();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFunc
pylith::materials::IsotropicLinearGenMaxwell::getKernelResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::f1v :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::f1v :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::f1v_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::f1v_refstate :
        NULL;

    PYLITH_METHOD_RETURN(f1u);
} // getKernelResidualStress


// ------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian G(t,s,\dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearGenMaxwell::getKernelJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::Jf3vu :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::Jf3vu :
        NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelJacobianElasticConstants


// ------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearGenMaxwell::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::cauchyStress :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::cauchyStress :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::cauchyStress_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::cauchyStress_refstate :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::IsotropicLinearGenMaxwell::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                    const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants"<<kernelConstants<<", dt="<<dt<<")");

    assert(kernelConstants);

    if (1 != kernelConstants->size()) { kernelConstants->resize(1);}
    (*kernelConstants)[0] = dt;

    PYLITH_METHOD_END;
} // updateKernelConstants


// ------------------------------------------------------------------------------------------------
// Add kernels for updating state variables.
void
pylith::materials::IsotropicLinearGenMaxwell::addKernelsUpdateStateVars(std::vector<ProjectKernels>* kernels,
                                                                        const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVars(kernels="<<kernels<<", coordsys="<<coordsys<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc funcViscousStrain =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::updateViscousStrain :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::updateViscousStrain :
        NULL;
    const PetscPointFunc funcTotalStrain =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::updateTotalStrain :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::updateTotalStrain :
        NULL;

    assert(kernels);
    size_t prevNumKernels = kernels->size();
    kernels->resize(prevNumKernels + 2);
    (*kernels)[prevNumKernels+0] = ProjectKernels("viscous_strain", funcViscousStrain);
    (*kernels)[prevNumKernels+1] = ProjectKernels("total_strain", funcTotalStrain);

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVars


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearGenMaxwell::getInterfaceKernelResidualF0Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelResidualF0Neg().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearGenMaxwell::getInterfaceKernelResidualF0Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelResidualF0Pos().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get f1 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearGenMaxwell::getInterfaceKernelResidualF1Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelResidualF1Neg().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get f1 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearGenMaxwell::getInterfaceKernelResidualF1Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelResidualF1Pos().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get Jf1lu kernel for LHS Jacobian F(t,s,dot{s}) for negative fault face.
PetscBdPointJac
pylith::materials::IsotropicLinearGenMaxwell::getInterfaceKernelJacobianF1Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelJacobianF1Neg().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get Jf1lu kernel for LHS Jacobian F(t,s,dot{s}) for positive fault face.
PetscBdPointJac
pylith::materials::IsotropicLinearGenMaxwell::getInterfaceKernelJacobianF1Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelJacobianF1Pos().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get Jf3lu kernel for LHS Jacobian F(t,s,dot{s}) for negative fault face.
PetscBdPointJac
pylith::materials::IsotropicLinearGenMaxwell::getInterfaceKernelJacobianF3Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelJacobianF3Neg().");
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get Jf3lu kernel for LHS Jacobian F(t,s,dot{s}) for positive fault face.
PetscBdPointJac
pylith::materials::IsotropicLinearGenMaxwell::getInterfaceKernelJacobianF3Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_COMPONENT_LOGICERROR(":TODO: @brad Implement getInterfaceKernelJacobianF3Pos().");
    return NULL;
}


// End of file
