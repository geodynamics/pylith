// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/materials/IsotropicLinearMaxwell.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryViscoelastic.hh" // USES AuxiliaryFactoryViscoelastic
#include "pylith/fekernels/IsotropicLinearMaxwell.hh" // USES IsotropicLinearMaxwell kernels
#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearMaxwell::IsotropicLinearMaxwell(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryViscoelastic),
    _useReferenceState(false) {
    _lhsJacobianTriggers = pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE;
    pylith::utils::PyreComponent::setName("isotropiclinearmaxwell");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearMaxwell::~IsotropicLinearMaxwell(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearMaxwell::deallocate(void) {
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearMaxwell::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearMaxwell::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicLinearMaxwell::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearMaxwell::addAuxiliarySubfields(void) {
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
    _auxiliaryFactory->addMaxwellTime();
    _auxiliaryFactory->addViscousStrain();
    _auxiliaryFactory->addTotalStrain();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFn*
pylith::materials::IsotropicLinearMaxwell::getKernelf1v(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1v(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* f1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f1v_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f1v_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f1v_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f1v_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(f1u);
} // getKernelf1v


// ------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJacFn*
pylith::materials::IsotropicLinearMaxwell::getKernelJf3vu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3vu(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::Jf3vu_infinitesimalStrain :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::Jf3vu_infinitesimalStrain :
        NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelJf3vu


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFn*
pylith::materials::IsotropicLinearMaxwell::getKernelf0Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0Neg(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFn* kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f0l_neg_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f0l_neg_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f0l_neg_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f0l_neg_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFn*
pylith::materials::IsotropicLinearMaxwell::getKernelf0Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFn* kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f0l_pos_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f0l_pos_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f0l_pos_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f0l_pos_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFn*
pylith::materials::IsotropicLinearMaxwell::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelCauchyStressVector(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::cauchyStress_infinitesimalStrain_asVector :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::cauchyStress_infinitesimalStrain_asVector :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::cauchyStress_infinitesimalStrain_refState_asVector :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::cauchyStress_infinitesimalStrain_refState_asVector :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelCauchyStressVector


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::IsotropicLinearMaxwell::updateKernelConstants(pylith::real_array* kernelConstants,
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
pylith::materials::IsotropicLinearMaxwell::addKernelsUpdateStateVars(std::vector<ProjectKernels>* kernels,
                                                                     const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVars(kernels="<<kernels<<", coordsys="<<coordsys<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFn* funcViscousStrain =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::viscousStrain_infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::viscousStrain_infinitesimalStrain_asVector :
        NULL;
    const PetscPointFn* funcTotalStrain =
        (3 == spaceDim) ? pylith::fekernels::Elasticity3D::infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector :
        NULL;

    assert(kernels);
    size_t prevNumKernels = kernels->size();
    kernels->resize(prevNumKernels + 2);
    (*kernels)[prevNumKernels+0] = ProjectKernels("viscous_strain", funcViscousStrain);
    (*kernels)[prevNumKernels+1] = ProjectKernels("total_strain", funcTotalStrain);

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVars


// End of file
