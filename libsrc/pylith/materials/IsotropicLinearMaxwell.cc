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

#include "pylith/materials/IsotropicLinearMaxwell.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryViscoelastic.hh" // USES AuxiliaryFactoryViscoelastic
#include "pylith/fekernels/IsotropicLinearMaxwell.hh" // USES IsotropicLinearMaxwell kernels
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearMaxwell::IsotropicLinearMaxwell(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryViscoelastic),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotropiclinearmaxwell");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearMaxwell::~IsotropicLinearMaxwell(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearMaxwell::deallocate(void) {
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearMaxwell::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearMaxwell::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicLinearMaxwell::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFunc
pylith::materials::IsotropicLinearMaxwell::getKernelResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f1v :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f1v :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f1v_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f1v_refstate :
        NULL;

    PYLITH_METHOD_RETURN(f1u);
} // getKernelResidualStress


// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearMaxwell::getKernelJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::Jf3vu :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::Jf3vu :
        NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelJacobianElasticConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearMaxwell::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::cauchyStress :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::cauchyStress :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::cauchyStress_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::cauchyStress_refstate :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables.
void
pylith::materials::IsotropicLinearMaxwell::addKernelsUpdateStateVars(std::vector<ProjectKernels>* kernels,
                                                                     const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVars(kernels="<<kernels<<", coordsys="<<coordsys<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc funcViscousStrain =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::updateViscousStrain :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::updateViscousStrain :
        NULL;
    const PetscPointFunc funcTotalStrain =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::updateTotalStrain :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::updateTotalStrain :
        NULL;

    assert(kernels);
    size_t prevNumKernels = kernels->size();
    kernels->resize(prevNumKernels + 2);
    (*kernels)[prevNumKernels+0] = ProjectKernels("viscous_strain", funcViscousStrain);
    (*kernels)[prevNumKernels+1] = ProjectKernels("total_strain", funcTotalStrain);

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVars


// End of file
