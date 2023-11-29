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

#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // implementation of object methods
#include "pylith/materials/AuxiliaryFactoryPoroelastic.hh" // USES AuxiliaryFactory

#include "pylith/fekernels/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearIncompElasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearPoroelasticity::IsotropicLinearPoroelasticity(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryPoroelastic),
    _useReferenceState(false),
    _useTensorPermeability(false) {
    pylith::utils::PyreComponent::setName("isotopiclinearporoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearPoroelasticity::~IsotropicLinearPoroelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearPoroelasticity::deallocate(void) {
    RheologyPoroelasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
void
pylith::materials::IsotropicLinearPoroelasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
bool
pylith::materials::IsotropicLinearPoroelasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use full tensor permeability?
void
pylith::materials::IsotropicLinearPoroelasticity::useTensorPermeability(const bool value) {
    PYLITH_COMPONENT_DEBUG("useTensorPermeability="<<value<<")");

    _useTensorPermeability = value;
} // useTensorPermeability


// ---------------------------------------------------------------------------------------------------------------------
// Use full tensor permeability?
bool
pylith::materials::IsotropicLinearPoroelasticity::useTensorPermeability(void) const {
    return _useTensorPermeability;
} // useTensorPermeability


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryPoroelastic*
pylith::materials::IsotropicLinearPoroelasticity::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearPoroelasticity::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress(); // numA - 7
        _auxiliaryFactory->addReferenceStrain(); // numA - 6
    } // if
    _auxiliaryFactory->addShearModulus(); // Shear Modulus, G, numA - 5
    _auxiliaryFactory->addDrainedBulkModulus(); // K_d, numA - 4
    _auxiliaryFactory->addBiotCoefficient(); // alpha, numA - 3
    _auxiliaryFactory->addBiotModulus(); // M, numA - 2
    if (_useTensorPermeability) {
        _auxiliaryFactory->addTensorPermeability(); // k, numA - 1
    } else {
        _auxiliaryFactory->addIsotropicPermeability(); // k, numA - 1
    }

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ================================ RHS ========================================

// ---------------------------------------------------------------------------------------------------------------------
// Select g0p function. Will only be used for the dynamic case.
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelg0p(const spatialdata::geocoords::CoordSys* coordsys,
                                                               const bool _useBodyForce,
                                                               const bool _gravityField,
                                                               const bool _useSourceDensity) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0p="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
    const int bitGravity = _gravityField ? 0x2 : 0x0;
    const int bitSourceDensity = _useSourceDensity ? 0x4 : 0x0;
    const int bitUse = bitBodyForce | bitGravity | bitSourceDensity;

    PetscPointFunc g0p = NULL;

    switch (bitUse) {
    case 0x0:
        g0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g0p :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p :
              NULL;
        break;
    case 0x1:
        g0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g0p :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p :
              NULL;
        break;
    case 0x2:
        g0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g0p :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p :
              NULL;
        break;
    case 0x4:
        g0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g0p_source :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_source :
              NULL; // aOff for sourceDensity is 3
        break;
    case 0x3:
        g0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g0p :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p :
              NULL;
        break;
    case 0x5:
        g0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g0p_source_body :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_source_body :
              NULL; // aOff for sourceDensity is 4
        break;
    case 0x6:
        g0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g0p_source_grav :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_source_grav :
              NULL; // aOff for sourceDensity is 4
        break;
    case 0x7:
        g0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g0p_source_grav_body :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_source_grav_body :
              NULL; // aOff for sourceDensity is 5
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for Poroelasticity RHS residual kernels.");
    } // switch

    PYLITH_METHOD_RETURN(g0p);
} // getKernelg0p


// ---------------------------------------------------------------------------------------------------------------------
// Get darcy velocity kernel
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelg1p_explicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                        const bool _gravityField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1p_implicit="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitTensorPermeability = _useTensorPermeability ? 0x1 : 0x0;
    const int bitGravityField = _gravityField ? 0x2 : 0x0;
    const int bitUse = bitTensorPermeability | bitGravityField;

    PetscPointFunc g1p = NULL;

    switch (bitUse) {
    case 0x0:
        g1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g1p :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p :
              NULL;
        break;
    case 0x1:
        g1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g1p_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_tensor_permeability :
              NULL;
        break;
    case 0x2:
        g1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g1p_gravity :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_gravity :
              NULL;
        break;
    case 0x3:
        g1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g1p_gravity_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_gravity_tensor_permeability :
              NULL;
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown combination of flags for  _useTensorPermeability="<<_useTensorPermeability<<", _gravityField="<<_gravityField<<").");
        throw std::logic_error("Unknown combination of flags.");
    } // switch

    PYLITH_METHOD_RETURN(g1p);
} // getKernelg1p_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1v_explicit(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitReferenceState = _useReferenceState ? 0x1 : 0x0;
    const int bitUse = bitReferenceState;

    PetscPointFunc g1v = NULL;

    switch (bitUse) {
    case 0x0:
        g1v = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g1v :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v :
              NULL;
        break;
    case 0x1:
        g1v = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g1v_refstate :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v_refstate :
              NULL;
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown combination of flags for getKernelg1v_explicit.");
        throw std::logic_error("Unknown combination of flags.");
    } // switch

    PYLITH_METHOD_RETURN(g1v);
} // getKernelResidualStress


// =============================== LHS =========================================

// ---------------------------------------------------------------------------------------------------------------------
// Get variation in fluid content kernel for LHS residual, F(t,s,\dot{s})
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelf0p_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0p_explicit="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f0p =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_explicit :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_explicit :
        NULL;

    PYLITH_METHOD_RETURN(f0p);
} // getKernelf0p_explicit


// ---------------------------------------------------------------------------------------------------------------------
// Select implicit f0p function.
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelf0p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                        const bool _useBodyForce,
                                                                        const bool _gravityField,
                                                                        const bool _useSourceDensity) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0p="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
    const int bitGravity = _gravityField ? 0x2 : 0x0;
    const int bitSourceDensity = _useSourceDensity ? 0x4 : 0x0;
    const int bitUse = bitBodyForce | bitGravity | bitSourceDensity;

    PetscPointFunc f0p = NULL;

    switch (bitUse) {
    case 0x0:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit :
              NULL;
        break;
    case 0x1:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit :
              NULL;
        break;
    case 0x2:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit :
              NULL;
        break;
    case 0x4:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit_source :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit_source :
              NULL; // aOff for sourceDensity is 3
        break;
    case 0x3:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit :
              NULL;
        break;
    case 0x5:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit_source_body :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit_source_body :
              NULL; // aOff for sourceDensity is 4
        break;
    case 0x6:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit_source_grav :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit_source_grav :
              NULL; // aOff for sourceDensity is 4
        break;
    case 0x7:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f0p_implicit_source_grav_body :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p_implicit_source_grav_body :
              NULL; // aOff for sourceDensity is 5
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for Poroelasticity LHS residual kernels.");
    } // switch

    PYLITH_METHOD_RETURN(f0p);
} // getKernelf0p_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1u(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitReferenceState = _useReferenceState ? 0x1 : 0x0;
    const int bitUse = bitReferenceState;

    PetscPointFunc f1u = NULL;

    switch (bitUse) {
    case 0x0:
        f1u = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1u :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1u :
              NULL;
        break;
    case 0x1:
        f1u = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1u_refstate :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1u_refstate :
              NULL;
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown combination of flags for getKernelf1u_implicit.");
        throw std::logic_error("Unknown combination of flags.");
    } // switch

    PYLITH_METHOD_RETURN(f1u);
} // getKernelf1u_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get darcy velocity kernel
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelf1p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                        const bool _useBodyForce,
                                                                        const bool _gravityField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1p_implicit="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitTensorPermeability = _useTensorPermeability ? 0x1 : 0x0;
    const int bitBodyForce = _useBodyForce ? 0x2 : 0x0;
    const int bitGravityField = _gravityField ? 0x4 : 0x0;
    const int bitUse = bitTensorPermeability | bitBodyForce | bitGravityField;

    PetscPointFunc f1p = NULL;

    switch (bitUse) {
    case 0x0:
        f1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1p :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1p :
              NULL;
        break;
    case 0x1:
        f1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1p_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1p_tensor_permeability :
              NULL;
        break;
    case 0x2:
        f1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1p_body :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1p_body :
              NULL;
        break;
    case 0x3:
        f1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1p_body_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1p_body_tensor_permeability :
              NULL;
        break;
    case 0x4:
        f1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1p_gravity :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1p_gravity :
              NULL;
        break;
    case 0x5:
        f1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1p_gravity_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1p_gravity_tensor_permeability :
              NULL;
        break;
    case 0x6:
        f1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1p_body_gravity :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1p_body_gravity :
              NULL;
        break;
    case 0x7:
        f1p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::f1p_body_gravity_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f1p_body_gravity_tensor_permeability :
              NULL;
        break;

    default:
        PYLITH_COMPONENT_ERROR("Unknown combination of flags for  _useTensorPermeability="<<_useTensorPermeability<<", _useBodyForce="<<_useBodyForce<<", _gravityField="<<_gravityField<<").");
        throw std::logic_error("Unknown combination of flags.");
    } // switch

    PYLITH_METHOD_RETURN(f1p);
} // getKernelf1p_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get poroelastic constants kernel for LHS Jacobian
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3uu(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf3uu :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf3uu :
        NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelLHSJacobianElasticConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for LHS Jacobian
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf2up(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf2up(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf2up =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf2up :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf2up :
        NULL;

    PYLITH_METHOD_RETURN(Jf2up);
} // getKernelJf2up


// ---------------------------------------------------------------------------------------------------------------------
// Get lambda kernel for LHS Jacobian
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf2ue(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf2ue(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf2ue =
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf2ue :
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf2ue :
        NULL;

    PYLITH_METHOD_RETURN(Jf2ue);
} // getKernelJf2ue


// ---------------------------------------------------------------------------------------------------------------------
// Get Specific storage kernel for LHS Jacobian F(t,s, \dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pp(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf0pp =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf0pp :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pp :
        NULL;

    PYLITH_METHOD_RETURN(Jf0pp);
} // getKernelJf0pp


// ---------------------------------------------------------------------------------------------------------------------
// Get Darcy Conductivity kernel for LHS Jacobian
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf3pp(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3pp(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3pp =
        (!_useTensorPermeability && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf3pp :
        (!_useTensorPermeability && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf3pp :
        (_useTensorPermeability && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf3pp_tensor_permeability :
        (_useTensorPermeability && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf3pp_tensor_permeability :
        NULL;

    PYLITH_METHOD_RETURN(Jf3pp);
} // getKerneJf3pp


// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for LHS Jacobian F(t,s, \dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pe(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf0pe =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf0pe :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pe :
        NULL;

    PYLITH_METHOD_RETURN(Jf0pe);
} // getKernelJf0pe


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf0ppdot(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0ppdot(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf0ppdot =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf0ppdot :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0ppdot :
        NULL;

    PYLITH_METHOD_RETURN(Jf0ppdot);
} // getKernelJf0ppdot


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf0pedot(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pedot(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf0pedot =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jf0pedot :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pedot :
        NULL;

    PYLITH_METHOD_RETURN(Jf0pedot);
} // getKernelJf0pedot


// =========================== DERIVED FIELDS ==================================

// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelCauchyStressVector(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ?  pylith::fekernels::IsotropicLinearPoroelasticity3D::cauchyStress_infinitesimalStrain_asVector :
        (!_useReferenceState && 2 == spaceDim) ?  pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::cauchyStress_infinitesimalStrain_asVector :
        (_useReferenceState && 3 == spaceDim) ?  pylith::fekernels::IsotropicLinearPoroelasticity3D::cauchyStress_infinitesimalStrain_refState_asVector :
        (_useReferenceState && 2 == spaceDim) ?  pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::cauchyStress_infinitesimalStrain_refState_asVector :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelCauchyStressVector

// ---------------------------------------------------------------------------------------------------------------------
// Get water content kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelWaterContent(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelWaterContent(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (3 == spaceDim) ?  pylith::fekernels::IsotropicLinearPoroelasticity3D::waterContent_asScalar :
        (2 == spaceDim) ?  pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::waterContent_asScalar :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelWaterContent

// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::IsotropicLinearPoroelasticity::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                        const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants"<<kernelConstants<<", dt="<<dt<<")");

    assert(kernelConstants);

    if (1 != kernelConstants->size()) { kernelConstants->resize(1);}
    (*kernelConstants)[0] = dt;

    PYLITH_METHOD_END;
} // updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables, implicit.
void
pylith::materials::IsotropicLinearPoroelasticity::addKernelsUpdateStateVarsImplicit(std::vector<ProjectKernels>* kernels,
                                                                                    const spatialdata::geocoords::CoordSys* coordsys,
                                                                                    const bool _useStateVars) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVarsImplicit(kernels="<<kernels<<", coordsys="<<coordsys<<")");
    if (_useStateVars) {
        const int spaceDim = coordsys->getSpaceDim();

        const PetscPointFunc funcPorosity =
            (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::updatePorosityImplicit :
            (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::updatePorosityImplicit :
            NULL;

        assert(kernels);
        size_t prevNumKernels = kernels->size();
        kernels->resize(prevNumKernels + 1);
        (*kernels)[prevNumKernels+0] = ProjectKernels("porosity", funcPorosity);
    }

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVarsImplicit


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables, explicit.
void
pylith::materials::IsotropicLinearPoroelasticity::addKernelsUpdateStateVarsExplicit(std::vector<ProjectKernels>* kernels,
                                                                                    const spatialdata::geocoords::CoordSys* coordsys,
                                                                                    const bool _useStateVars) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVarsExplicit(kernels="<<kernels<<", coordsys="<<coordsys<<")");
    if (_useStateVars) {
        const int spaceDim = coordsys->getSpaceDim();

        const PetscPointFunc funcPorosity =
            (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::updatePorosityExplicit :
            (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::updatePorosityExplicit :
            NULL;

        assert(kernels);
        size_t prevNumKernels = kernels->size();
        kernels->resize(prevNumKernels + 1);
        (*kernels)[prevNumKernels+0] = ProjectKernels("porosity", funcPorosity);
    }

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVarsExplicit


// End of file
