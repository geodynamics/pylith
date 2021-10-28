// -*- C++ -*-
//
// ------------------------------------------------------------------------------------------------
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
// ------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "Material.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(void) :
    _gravityField(NULL),
    _materialId(0),
    _descriptiveLabel("") {}


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
// Set value of label material-id used to identify material cells.
void
pylith::materials::Material::setMaterialId(const int value) {
    PYLITH_COMPONENT_DEBUG("setMaterialId(value="<<value<<")");

    _materialId = value;
} // setMaterialId


// ------------------------------------------------------------------------------------------------
// Get value of label material-id used to identify material cells.
int
pylith::materials::Material::getMaterialId(void) const {
    return _materialId;
} // getMaterialId


// ------------------------------------------------------------------------------------------------
// Set descriptive label of material.
void
pylith::materials::Material::setDescriptiveLabel(const char* value) {
    PYLITH_COMPONENT_DEBUG("setDescriptiveLabel(value="<<value<<")");

    _descriptiveLabel = value;
} // setDescriptiveLabel


// ------------------------------------------------------------------------------------------------
// Get label of material.
const char*
pylith::materials::Material::getDescriptiveLabel(void) const {
    return _descriptiveLabel.c_str();
} // getDescriptiveLabel


// ------------------------------------------------------------------------------------------------
// Set gravity field.
void
pylith::materials::Material::setGravityField(spatialdata::spatialdb::GravityField* const g) {
    _gravityField = g;
} // setGravityField


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
