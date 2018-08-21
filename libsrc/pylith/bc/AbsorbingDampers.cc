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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AbsorbingDampers.hh"// implementation of object methods

#include "AbsorbingDampersAuxiliaryFactory.hh"// USES AuxiliaryFactory

#include "pylith/fekernels/AbsorbingDampers.hh"// USES AbsorbingDampers kernels

#include "pylith/feassemble/IntegratorBoundary.hh"// USES IntegratorBoundary
#include "pylith/topology/Mesh.hh"// USES Mesh
#include "pylith/topology/Field.hh"// USES Field
#include "pylith/topology/FieldOps.hh"// USES FieldOps
#include "spatialdata/units/Nondimensional.hh"// USES Nondimensional

#include "pylith/utils/journals.hh"// USES PYLITH_COMPONENT_*

#include <strings.h>// USES strcasecmp()
#include <cassert>// USES assert()
#include <stdexcept>// USES std::runtime_error
#include <sstream>// USES std::ostringstream
#include <typeinfo>// USES typeid

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorBoundary::ResidualKernels ResidualKernels;

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace bc {
        class _AbsorbingDampers {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /** Set kernels for RHS residual.
             *
             * @param[out] integrator Integrator for boundary condition.
             * @param[in] solution Solution field.
             */
            static
            void setKernelsRHSResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                       const pylith::topology::Field& solution);

            static const char* subfieldName;
        };

        // _NeumannTimeDependent
        const char* _AbsorbingDampers::subfieldName = "velocity";

    }// bc
}// pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampers::AbsorbingDampers(void) :
    _auxiliaryFactory(new pylith::bc::AbsorbingDampersAuxiliaryFactory),
    _boundaryLabel("") {
    PyreComponent::name("absorbingdampers");

    _refDir1[0] = 0.0;
    _refDir1[1] = 0.0;
    _refDir1[2] = 1.0;

    _refDir2[0] = 0.0;
    _refDir2[1] = 1.0;
    _refDir2[2] = 0.0;
}// constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampers::~AbsorbingDampers(void) {
    deallocate();
}// destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::AbsorbingDampers::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    Physics::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;

    PYLITH_METHOD_END;
}// deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set label marking boundary associated with boundary condition surface.
void
pylith::bc::AbsorbingDampers::setMarkerLabel(const char* value) {
    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition label.");
    }// if

    _boundaryLabel = value;
}// setMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label marking boundary associated with boundary condition surface.
const char*
pylith::bc::AbsorbingDampers::getMarkerLabel(void) const {
    return _boundaryLabel.c_str();
}// getMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Set first choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::bc::AbsorbingDampers::setRefDir1(const PylithReal vec[3]) {
    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 1 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for boundary condition '" << _boundaryLabel << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    }// if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    }// for
}// setRefDir1


// ---------------------------------------------------------------------------------------------------------------------
// Set second choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::bc::AbsorbingDampers::setRefDir2(const PylithReal vec[3]) {
    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 2 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for boundary condition '" << _boundaryLabel << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    }// if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    }// for
}// setRefDir2


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    if (!solution.hasSubfield(_AbsorbingDampers::subfieldName)) {
        std::ostringstream msg;
        msg << "Cannot apply Neumann boundary condition to field '"<< _AbsorbingDampers::subfieldName
            << "'; field is not in solution.";
        throw std::runtime_error(msg.str());
    }// if

    const pylith::topology::Field::SubfieldInfo& info = solution.subfieldInfo(_AbsorbingDampers::subfieldName);
    if (pylith::topology::Field::VECTOR != info.description.vectorFieldType) {
        std::ostringstream msg;
        msg << "Absorbing boundary condition cannot be applied to non-vector field '"<< _AbsorbingDampers::subfieldName << "' in solution.";
        throw std::runtime_error(msg.str());
    }// if

    PYLITH_METHOD_END;
}// verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::bc::AbsorbingDampers::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.label()<<")");

    pylith::feassemble::IntegratorBoundary* integrator = new pylith::feassemble::IntegratorBoundary(this);assert(integrator);
    integrator->setMarkerLabel(getMarkerLabel());

    _AbsorbingDampers::setKernelsRHSResidual(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
}// createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::bc::AbsorbingDampers::createConstraint(const pylith::topology::Field& solution) {
    return NULL;
}// createConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::bc::AbsorbingDampers::createAuxiliaryField(const pylith::topology::Field& solution,
                                                   const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.label()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->label("AbsorbingDampers auxiliary field");

    assert(_auxiliaryFactory);
    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.dimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.

    _auxiliaryFactory->addDensity();// 0
    _auxiliaryFactory->addVp();// 1
    _auxiliaryFactory->addVs();// 2

    auxiliaryField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->zeroLocal();

    assert(_auxiliaryFactory);
    _auxiliaryFactory->initializeSubfields();

    PYLITH_METHOD_RETURN(auxiliaryField);
}// createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::bc::AbsorbingDampers::createDerivedField(const pylith::topology::Field& solution,
                                                 const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.label()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
}// createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::AbsorbingDampers::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
}// _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::bc::AbsorbingDampers::_updateKernelConstants(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelConstants(dt="<<dt<<")");

    if (6 != _kernelConstants.size()) { _kernelConstants.resize(6);}
    _kernelConstants[0] = _refDir1[0];
    _kernelConstants[1] = _refDir1[1];
    _kernelConstants[2] = _refDir1[2];
    _kernelConstants[3] = _refDir2[0];
    _kernelConstants[4] = _refDir2[1];
    _kernelConstants[5] = _refDir2[2];

    PYLITH_METHOD_END;
}// _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::bc::_AbsorbingDampers::setKernelsRHSResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                                     const topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    //PYLITH_COMPONENT_DEBUG("_setKernelsRHSResidual(integrator="<<integrator<<", solution="<<solution.label()<<")");

    PetscBdPointFunc g0 = pylith::fekernels::AbsorbingDampers::g0;
    PetscBdPointFunc g1 = NULL;

    std::vector<ResidualKernels> kernels(1);
    kernels[0] = ResidualKernels(subfieldName, g0, g1);

    assert(integrator);
    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END;
}// _setKernelsRHSResidual


// End of file
