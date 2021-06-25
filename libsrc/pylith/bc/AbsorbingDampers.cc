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

#include "AbsorbingDampers.hh" // implementation of object methods

#include "AbsorbingDampersAuxiliaryFactory.hh" // USES AuxiliaryFactory

#include "pylith/fekernels/AbsorbingDampers.hh" // USES AbsorbingDampers kernels

#include "pylith/feassemble/IntegratorBoundary.hh" // USES IntegratorBoundary
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES typeid

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
             * @param[in] bc Absorbing dampers boundary condition.
             * @param[in] solution Solution field.
             */
            static
            void setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                    const pylith::bc::AbsorbingDampers& bc,
                                    const pylith::topology::Field& solution);

            static const char* pyreComponent;
        };
        const char* _AbsorbingDampers::pyreComponent = "absorbingdampers";

    } // bc
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::AbsorbingDampers::AbsorbingDampers(void) :
    _auxiliaryFactory(new pylith::bc::AbsorbingDampersAuxiliaryFactory),
    _boundaryLabel("") {
    PyreComponent::setName(_AbsorbingDampers::pyreComponent);

    _subfieldName = "velocity";
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::AbsorbingDampers::~AbsorbingDampers(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::AbsorbingDampers::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    BoundaryCondition::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::AbsorbingDampers::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    BoundaryCondition::verifyConfiguration(solution);

    const pylith::topology::Field::SubfieldInfo& info = solution.getSubfieldInfo(_subfieldName.c_str());
    if (pylith::topology::Field::VECTOR != info.description.vectorFieldType) {
        std::ostringstream msg;
        msg << "Absorbing boundary condition cannot be applied to non-vector field '"<< _subfieldName << "' in solution.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::bc::AbsorbingDampers::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorBoundary* integrator = new pylith::feassemble::IntegratorBoundary(this);assert(integrator);
    integrator->setMarkerLabel(getMarkerLabel());
    integrator->setSubfieldName(getSubfieldName());
    integrator->setLabelName(getMarkerLabel());

    _AbsorbingDampers::setKernelsResidual(integrator, *this, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::bc::AbsorbingDampers::createConstraint(const pylith::topology::Field& solution) {
    return NULL;
} // createConstraint


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::bc::AbsorbingDampers::createAuxiliaryField(const pylith::topology::Field& solution,
                                                   const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("AbsorbingDampers auxiliary field");

    assert(_auxiliaryFactory);
    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.getDimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.

    _auxiliaryFactory->addDensity(); // 0
    _auxiliaryFactory->addVp(); // 1
    _auxiliaryFactory->addVs(); // 2

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    assert(_auxiliaryFactory);
    _auxiliaryFactory->setValuesFromDB();

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::bc::AbsorbingDampers::createDerivedField(const pylith::topology::Field& solution,
                                                 const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary field factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::AbsorbingDampers::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


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
} // _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::bc::_AbsorbingDampers::setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                                  const pylith::bc::AbsorbingDampers& bc,
                                                  const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_AbsorbingDampers::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "_AbsorbingDampers::_setKernelsRHSResidual(integrator="<<integrator<<", bc="<<typeid(bc).name()
          <<", solution="<<solution.getLabel()<<")"<<pythia::journal::endl;

    PetscBdPointFunc g0 = pylith::fekernels::AbsorbingDampers::g0;
    PetscBdPointFunc g1 = NULL;

    std::vector<ResidualKernels> kernels(1);
    kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::RESIDUAL_RHS, g0, g1);

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsRHSResidual


// End of file
