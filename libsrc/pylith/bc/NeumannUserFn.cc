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

#include "NeumannUserFn.hh" // implementation of object methods

#include "pylith/fekernels/NeumannTimeDependent.hh" // USES NeumannTimeDepndent kernels

#include "pylith/feassemble/IntegratorBoundary.hh" // USES IntegratorBoundary
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorBoundary::ResidualKernels ResidualKernels;

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace bc {
        class _NeumannUserFn {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            /** Set kernels for RHS residual.
             *
             * @param[out] integrator Integrator for boundary condition.
             * @param[in] bc Neumann time-dependent boundary condition.
             * @param[in] solution Solution field.
             * @param[in] formulation Formulation for equations.
             */
            static
            void setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                    const pylith::bc::NeumannUserFn& bc,
                                    const pylith::topology::Field& solution,
                                    const pylith::problems::Physics::FormulationEnum formulation);

            static const char* pyreComponent;

        }; // _NeumannUserFn
        const char* _NeumannUserFn::pyreComponent = "neumannuserfn";

    } // bc
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::NeumannUserFn::NeumannUserFn(void) :
    _fn(NULL) {
    PyreComponent::setName("neumannuserfn");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::NeumannUserFn::~NeumannUserFn(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::NeumannUserFn::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    BoundaryCondition::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set user function specifying field on boundary.
void
pylith::bc::NeumannUserFn::setUserFn(PetscBdPointFunc fn) {
    PYLITH_COMPONENT_DEBUG("setUserFn(fn="<<fn<<")");

    _fn = fn;
} // setUserFn


// ------------------------------------------------------------------------------------------------
// Get user function specifying field on boundary.
PetscBdPointFunc
pylith::bc::NeumannUserFn::getUserFn(void) const {
    return _fn;
} // getUserFn


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::bc::NeumannUserFn::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorBoundary* integrator = new pylith::feassemble::IntegratorBoundary(this);assert(integrator);
    integrator->setMarkerLabel(getMarkerLabel());
    integrator->setSubfieldName(getSubfieldName());
    integrator->setLabelName(getMarkerLabel());

    _NeumannUserFn::setKernelsResidual(integrator, *this, solution, _formulation);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<pylith::feassemble::Constraint*>
pylith::bc::NeumannUserFn::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getLabel()<<") empty method");
    std::vector<pylith::feassemble::Constraint*> constraintArray;

    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::bc::NeumannUserFn::createAuxiliaryField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("NeumannUserFn auxiliary field (not used)");

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying auxiliary field");
        auxiliaryField->view("Neumann auxiliary field");
    } // if

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::bc::NeumannUserFn::createDerivedField(const pylith::topology::Field& solution,
                                              const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::bc::NeumannUserFn::_getAuxiliaryFactory(void) {
    return NULL;
} // _getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::bc::NeumannUserFn::_updateKernelConstants(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateKernelConstants(dt="<<dt<<")");

    if (6 != _kernelConstants.size()) { _kernelConstants.resize(6);}
    _kernelConstants[0] = _refDir1[0];
    _kernelConstants[1] = _refDir1[1];
    _kernelConstants[2] = _refDir1[2];
    _kernelConstants[3] = _refDir2[0];
    _kernelConstants[4] = _refDir2[1];
    _kernelConstants[5] = _refDir2[2];

    PYLITH_METHOD_END;
} // _updateKernelConstants


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::bc::_NeumannUserFn::setKernelsResidual(pylith::feassemble::IntegratorBoundary* integrator,
                                               const pylith::bc::NeumannUserFn& bc,
                                               const topology::Field& solution,
                                               const pylith::problems::Physics::FormulationEnum formulation) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug(_NeumannUserFn::pyreComponent);
    debug << pythia::journal::at(__HERE__)
          << "setKernelsResidual(integrator="<<integrator<<", bc="<<typeid(bc).name()<<", solution="
          << solution.getLabel()<<")"
          << pythia::journal::endl;

    PetscBdPointFunc r0 = bc.getUserFn();
    PetscBdPointFunc r1 = NULL;

    std::vector<ResidualKernels> kernels(1);
    switch (formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::RESIDUAL_LHS, r0, r1);
        break;
    case pylith::problems::Physics::DYNAMIC_IMEX:
    case pylith::problems::Physics::DYNAMIC:
        kernels[0] = ResidualKernels(bc.getSubfieldName(), pylith::feassemble::Integrator::RESIDUAL_RHS, r0, r1);
        break;
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown formulation for equations ("<<formulation<<").");
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // setKernelsResidual


// End of file
