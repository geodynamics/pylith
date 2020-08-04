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

#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // implementation of object methods
#include "pylith/materials/AuxiliaryFactoryPoroelastic.hh" // USES AuxiliaryFactory

#include "pylith/fekernels/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearIncompElasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo>  // USES typeid()

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
    _auxiliaryFactory->addUndrainedBulkModulus(); // K_u, numA - 4
    _auxiliaryFactory->addBiotCoefficient();  // alpha, numA - 3
    _auxiliaryFactory->addBiotModulus(); // M, numA - 2
    if (_useTensorPermeability) {
        _auxiliaryFactory->addTensorPermeability(); // k, numA - 1
    } else {
        _auxiliaryFactory->addIsotropicPermeability();  // k, numA - 1
    }


    PYLITH_METHOD_END;
} // addAuxiliarySubfields

// ================================ RHS ========================================
// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelg1u(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1u(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointFunc g1u =   (!_useReferenceState) ?
                           pylith::fekernels::IsotropicLinearPoroelasticity::g1u :
                           pylith::fekernels::IsotropicLinearPoroelasticity::g1u_refstate;


    PYLITH_METHOD_RETURN(g1u);
} // getKernelg1u

// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelg1v(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1v(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointFunc g1v =   (!_useReferenceState) ?
                           pylith::fekernels::IsotropicLinearPoroelasticity::g1v :
                           pylith::fekernels::IsotropicLinearPoroelasticity::g1v_refstate;

    PYLITH_METHOD_RETURN(g1v);
} // getKernelg1v


// ---------------------------------------------------------------------------------------------------------------------
// Get darcy velocity kernel for RHS residual, G(t,s)
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelg1p(const spatialdata::geocoords::CoordSys* coordsys, const bool _gravityField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1p="<<typeid(coordsys).name()<<")");

    PetscPointFunc g1p =
        (!_gravityField && !_useTensorPermeability) ? pylith::fekernels::IsotropicLinearPoroelasticity::g1p_NoGrav :
        (!_gravityField && _useTensorPermeability) ? pylith::fekernels::IsotropicLinearPoroelasticity::g1p_NoGrav_tensorPerm :
        (_gravityField && !_useTensorPermeability) ? pylith::fekernels::IsotropicLinearPoroelasticity::g1p_Grav :
        (_gravityField && _useTensorPermeability) ?  pylith::fekernels::IsotropicLinearPoroelasticity::g1p_NoGrav_tensorPerm :
        NULL;

    PYLITH_METHOD_RETURN(g1p);
  } // getKernelg1p


// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJg3uu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJg3uu(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg3uu = pylith::fekernels::IsotropicLinearPoroelasticity::Jg3uu;

    PYLITH_METHOD_RETURN(Jg3uu);
} // getKernelRHSJacobianElasticConstants

// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJg3vu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJg3vu(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg3vu = pylith::fekernels::IsotropicLinearPoroelasticity::Jg3vu;

    PYLITH_METHOD_RETURN(Jg3vu);
} // getKernelJg3vu

// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJg2up(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJg2up(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg2up = pylith::fekernels::IsotropicLinearPoroelasticity::Jg2up;

    PYLITH_METHOD_RETURN(Jg2up);
} // getKernelJg2up

// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJg2vp(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSJacobianBiotCoefficient(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg2vp = pylith::fekernels::IsotropicLinearPoroelasticity::Jg2vp;

    PYLITH_METHOD_RETURN(Jg2vp);
} // getKerneJg2vp

// ---------------------------------------------------------------------------------------------------------------------
// Get lambda kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJg2ue(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJg2ue(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg2ue = pylith::fekernels::IsotropicLinearPoroelasticity::Jg2ue;

    PYLITH_METHOD_RETURN(Jg2ue);
} // getKernelJg2ue


// ---------------------------------------------------------------------------------------------------------------------
// Get Darcy Conductivity kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJg3pp(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJg3pp(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg3pp =
        (!_useTensorPermeability) ? pylith::fekernels::IsotropicLinearPoroelasticity::Jg3pp :
        (_useTensorPermeability) ? pylith::fekernels::IsotropicLinearPoroelasticity::Jg3pp_tensorPerm :
        NULL;

    PYLITH_METHOD_RETURN(Jg3pp);
} // getKerneJg3pp

// =============================== LHS =========================================

// ---------------------------------------------------------------------------------------------------------------------
// Get variation in fluid content kernel for LHS residual, F(t,s,\dot{s})
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelf0p(const spatialdata::geocoords::CoordSys* coordsys, const bool _useInertia) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelLHSVariationInFluidContent="<<typeid(coordsys).name()<<")");

    PetscPointFunc f0p = (!_useInertia) ?
                          pylith::fekernels::IsotropicLinearPoroelasticity::f0p_QS :
                          pylith::fekernels::IsotropicLinearPoroelasticity::f0p_DYN;

    PYLITH_METHOD_RETURN(f0p);
  } // getKernelf0p

// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for LHS Jacobian F(t,s, \dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const {
  PYLITH_METHOD_BEGIN;
  PYLITH_COMPONENT_DEBUG("getKernelJf0pe(coordsys="<<typeid(coordsys).name()<<")");

  PetscPointJac Jf0pe = pylith::fekernels::IsotropicLinearPoroelasticity::Jf0pe;

  PYLITH_METHOD_RETURN(Jf0pe);
} // getKernelJf0pe


// ---------------------------------------------------------------------------------------------------------------------
// Get Specific storage kernel for LHS Jacobian F(t,s, \dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const {
  PYLITH_METHOD_BEGIN;
  PYLITH_COMPONENT_DEBUG("getKernelJf0pp(coordsys="<<typeid(coordsys).name()<<")");

  PetscPointJac Jf0pp = pylith::fekernels::IsotropicLinearPoroelasticity::Jf0pp;

  PYLITH_METHOD_RETURN(Jf0pp);
} // getKernelJf0pp


// =========================== DERIVED FIELDS ==================================

// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel = (!_useReferenceState) ?
                              pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress :
                              pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refstate;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// End of file
