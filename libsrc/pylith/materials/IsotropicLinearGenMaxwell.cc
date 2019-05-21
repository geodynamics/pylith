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

#include "pylith/materials/IsotropicLinearGenMaxwell.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryViscoelastic.hh" // USES AuxiliaryFactoryViscoelastic
#include "pylith/fekernels/IsotropicLinearGenMaxwell.hh" // USES IsotropicLinearGenMaxwell kernels
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearGenMaxwell::IsotropicLinearGenMaxwell(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryViscoelastic),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotropiclineargenmaxwell");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearGenMaxwell::~IsotropicLinearGenMaxwell(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearGenMaxwell::deallocate(void) {
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearGenMaxwell::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearGenMaxwell::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicLinearGenMaxwell::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
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
    _auxiliaryFactory->addViscousStrainGeneralizedMaxwell(); // 4
    _auxiliaryFactory->addTotalStrain();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearGenMaxwell::getKernelRHSResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointFunc g1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g1v :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::g1v :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::g1v_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::g1v_refstate :
        NULL;

    PYLITH_METHOD_RETURN(g1u);
} // getKernelRHSResidualStress


// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearGenMaxwell::getKernelRHSJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointJac Jg3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::Jg3vu :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::Jg3vu :
        NULL;

    PYLITH_METHOD_RETURN(Jg3uu);
} // getKernelRHSJacobianElasticConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearGenMaxwell::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::stress :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::stress :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwell3D::stress_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain::stress_refstate :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables.
void
pylith::materials::IsotropicLinearGenMaxwell::addKernelsUpdateStateVars(std::vector<ProjectKernels>* kernels,
                                                                        const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVars(kernels="<<kernels<<", coordsys="<<coordsys<<")");

    const int spaceDim = coordsys->spaceDim();
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


// End of file
