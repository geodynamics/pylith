// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
//

#include <portinfo>

#include "BoundaryCondition.hh" // implementation of object methods

#include "pylith/bc/DiagnosticFieldFactory.hh" // USES DiagnosticFieldFactory
#include "pylith/feassemble/IntegratorBoundary.hh" // USES IntegratorBoundary
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/fekernels/BoundaryDirections.hh" // USES BoundaryDirections

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::createOutputLabel()

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cstring> // USES strlen()
#include <stdexcept> // USES std::runtime_error()

// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::BoundaryCondition::BoundaryCondition(void) :
    _subfieldName(""),
    _diagnosticFactory(new pylith::bc::DiagnosticFieldFactory) {
    _refDir1[0] = 0.0;
    _refDir1[1] = 0.0;
    _refDir1[2] = 1.0;

    _refDir2[0] = 0.0;
    _refDir2[1] = 1.0;
    _refDir2[2] = 0.0;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::BoundaryCondition::~BoundaryCondition(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::BoundaryCondition::deallocate(void) {
    Physics::deallocate();

    delete _diagnosticFactory;_diagnosticFactory = NULL;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of solution subfield associated with boundary condition.
void
pylith::bc::BoundaryCondition::setSubfieldName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setSubfieldName(value="<<value<<")");

    if (!value || (0 == strlen(value))) {
        std::ostringstream msg;
        msg << "Empty string given for name of solution subfield for boundary condition '"
            << getLabelName() <<"'.";
        throw std::runtime_error(msg.str());
    } // if
    _subfieldName = value;
} // setSubfieldName


// ------------------------------------------------------------------------------------------------
// Get name of solution subfield associated with boundary condition.
const char*
pylith::bc::BoundaryCondition::getSubfieldName(void) const {
    return _subfieldName.c_str();
} // getSubfieldName


// ------------------------------------------------------------------------------------------------
// Set first choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::bc::BoundaryCondition::setRefDir1(const PylithReal vec[3]) {
    PYLITH_COMPONENT_DEBUG("setRefDir1(vec="<<vec[0]<<","<<vec[1]<<","<<vec[2]<<")");

    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 1 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for boundary condition '" << getLabelName() << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    } // for
} // setRefDir1


// ------------------------------------------------------------------------------------------------
// Set second choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::bc::BoundaryCondition::setRefDir2(const PylithReal vec[3]) {
    PYLITH_COMPONENT_DEBUG("setRefDir2(vec="<<vec[0]<<","<<vec[1]<<","<<vec[2]<<")");

    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 2 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]
            <<") for boundary condition '" << getLabelName() << "' is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir2[i] = vec[i] / mag;
    } // for
} // setRefDir2


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::BoundaryCondition::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    if (!solution.hasSubfield(_subfieldName.c_str())) {
        std::ostringstream msg;
        msg << "Cannot apply boundary condition to field '"<< _subfieldName
            << "'; field is not in solution.";
        throw std::runtime_error(msg.str());
    } // if

    const PetscDM dmSoln = solution.getDM();
    PetscBool hasLabel = PETSC_FALSE;
    PetscErrorCode err = DMHasLabel(dmSoln, getLabelName(), &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << getLabelName() << "' in boundary condition '"
            << PyreComponent::getIdentifier() <<"'.";
        throw std::runtime_error(msg.str());
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create diagnostic field.
pylith::topology::Field*
pylith::bc::BoundaryCondition::createDiagnosticField(const pylith::topology::Field& solution,
                                                     const pylith::topology::Mesh& physicsMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDiagnosticField(solution="<<solution.getLabel()<<", physicsMesh=)"<<typeid(physicsMesh).name()<<")");

    assert(_normalizer);

    pylith::topology::Field* diagnosticField = new pylith::topology::Field(physicsMesh);assert(diagnosticField);
    diagnosticField->setLabel("BoundaryCondition diagnostic field");
    pylith::topology::FieldOps::createOutputLabel(diagnosticField);

    assert(_diagnosticFactory);
    const pylith::topology::FieldBase::Discretization& discretization = solution.getSubfieldInfo("displacement").fe;
    const PylithInt cellDim = solution.getSpaceDim()-1;
    const bool isFaultOnly = false;
    _diagnosticFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder, cellDim,
                                                  isFaultOnly, discretization.cellBasis, discretization.feSpace,
                                                  discretization.isBasisContinuous);

    assert(_diagnosticFactory);
    assert(_normalizer);
    _diagnosticFactory->initialize(diagnosticField, *_normalizer, solution.getSpaceDim());

    _diagnosticFactory->addNormalDir(); // 0
    _diagnosticFactory->addTangentialDirHoriz(); // 1
    if (solution.getSpaceDim() > 2) {
        _diagnosticFactory->addTangentialDirVert(); // 2
    } // if

    diagnosticField->subfieldsSetup();
    diagnosticField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *diagnosticField);
    diagnosticField->allocate();
    diagnosticField->createOutputVector();

    PYLITH_METHOD_RETURN(diagnosticField);
}


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::bc::BoundaryCondition::_updateKernelConstants(const PylithReal dt) {
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


// ------------------------------------------------------------------------------------------------
// Set kernels for computing diagnostic field.
void
pylith::bc::BoundaryCondition::_setKernelsDiagnosticField(pylith::feassemble::IntegratorBoundary* integrator,
                                                          const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDiagnosticField(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");
    typedef pylith::feassemble::IntegratorBoundary::ProjectKernels ProjectKernels;

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    const PylithInt spaceDim = solution.getSpaceDim();
    std::vector<ProjectKernels> kernels(spaceDim);
    kernels[0] = ProjectKernels("normal_dir", pylith::fekernels::BoundaryDirections::normalDir);
    if (spaceDim == 2) {
        kernels[1] = ProjectKernels("tangential_dir", pylith::fekernels::BoundaryDirections::tangentialDirHoriz);
    } else {
        kernels[1] = ProjectKernels("horizontal_tangential_dir", pylith::fekernels::BoundaryDirections::tangentialDirHoriz);
        kernels[2] = ProjectKernels("vertical_tangential_dir", pylith::fekernels::BoundaryDirections::tangentialDirVert);
    } // if/else

    assert(integrator);
    integrator->setKernelsDiagnosticField(kernels);

    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Set kernels for computing diagnostic field.
void
pylith::bc::BoundaryCondition::_setKernelsDiagnosticField(pylith::feassemble::Constraint* constraint,
                                                          const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDiagnosticField(constraint="<<constraint<<", solution="<<solution.getLabel()<<")");
    typedef pylith::feassemble::Constraint::ProjectKernels ProjectKernels;

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    const PylithInt spaceDim = solution.getSpaceDim();
    std::vector<ProjectKernels> kernels(spaceDim);
    kernels[0] = ProjectKernels("normal_dir", pylith::fekernels::BoundaryDirections::normalDir);
    if (spaceDim == 2) {
        kernels[1] = ProjectKernels("tangential_dir", pylith::fekernels::BoundaryDirections::tangentialDirHoriz);
    } else {
        kernels[1] = ProjectKernels("horizontal_tangential_dir", pylith::fekernels::BoundaryDirections::tangentialDirHoriz);
        kernels[2] = ProjectKernels("vertical_tangential_dir", pylith::fekernels::BoundaryDirections::tangentialDirVert);
    } // if/else

    assert(constraint);
    constraint->setKernelsDiagnosticField(kernels);

    PYLITH_METHOD_END;
}


// End of file
