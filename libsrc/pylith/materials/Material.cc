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

#include "pylith/materials/Material.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(void) :
    _gravityField(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::Material::~Material(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Material::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::problems::Physics::deallocate();

    _gravityField = NULL; // :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set gravity field.
void
pylith::materials::Material::setGravityField(spatialdata::spatialdb::GravityField* const g) {
    _gravityField = g;
} // setGravityField


// ------------------------------------------------------------------------------------------------
// Set body force kernels (intended for use in MMS tests).
void
pylith::materials::Material::setMMSBodyForceKernels(const std::vector<pylith::feassemble::IntegratorDomain::ResidualKernels>& kernels) {
    _mmsBodyForceKernels = kernels;
}


// ------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<pylith::feassemble::Constraint*>
pylith::materials::Material::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getLabel()<<") empty method");
    std::vector<pylith::feassemble::Constraint*> constraintArray;

    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Get default PETSc solver options appropriate for material.
pylith::utils::PetscOptions*
pylith::materials::Material::getSolverDefaults(const bool isParallel,
                                               const bool hasFault) const {
    return NULL;
}


// ------------------------------------------------------------------------------------------------
// Get residual kernels for an interior interface bounding material.
std::vector<pylith::materials::Material::InterfaceResidualKernels>
pylith::materials::Material::getInterfaceKernelsResidual(const pylith::topology::Field& solution,
                                                         pylith::feassemble::IntegratorInterface::FaceEnum face) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceResidualKernels(solution="<<solution.getLabel()<<", face="<<face<<") empty method");
    std::vector<InterfaceResidualKernels> kernels;

    PYLITH_METHOD_RETURN(kernels);
} // getInterfaceKernelsResidual


// ------------------------------------------------------------------------------------------------
// Get Jacobian kernels for an interior interface bounding material.
std::vector<pylith::materials::Material::InterfaceJacobianKernels>
pylith::materials::Material::getInterfaceKernelsJacobian(const pylith::topology::Field& solution,
                                                         pylith::feassemble::IntegratorInterface::FaceEnum face) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceJacobianKernels(solution="<<solution.getLabel()<<", face="<<face<<") empty method");
    std::vector<InterfaceJacobianKernels> kernels;

    PYLITH_METHOD_RETURN(kernels);
} // getInterfaceKernelsJacobian


// End of file
