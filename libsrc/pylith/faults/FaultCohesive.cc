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

#include "pylith/faults/FaultCohesive.hh" // implementation of object methods

#include "pylith/faults/TopologyOps.hh" // USES TopologyOps
#include "pylith/faults/AuxiliaryFieldFactory.hh" // USES AuxiliaryFieldFactory
#include "pylith/faults/DiagnosticFieldFactory.hh" // USES DiagnosticFieldFactory
#include "pylith/faults/DerivedFieldFactory.hh" // USES DerivedFieldFactory
#include "pylith/fekernels/BoundaryDirections.hh" // USES BoundaryDirections
#include "pylith/fekernels/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/InterfacePatches.hh" // USES InterfacePatches
#include "pylith/feassemble/ConstraintSimple.hh" // USES ConstraintSimple
#include "pylith/feassemble/FEKernelKey.hh" // USES FEKernelKey
#include "pylith/topology/Mesh.hh" // USES Mesh::cells_label_name
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDiscretization()
#include "pylith/topology/MeshOps.hh" // USES MeshOps::checkTopology()
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <utility> // USES std::pair
#include <map> // USES std::map
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorInterface::ProjectKernels ProjectKernels;

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesive::FaultCohesive(void) :
    _auxiliaryFactory(new pylith::faults::AuxiliaryFieldFactory),
    _diagnosticFactory(new pylith::faults::DiagnosticFieldFactory),
    _derivedFactory(new pylith::faults::DerivedFieldFactory),
    _surfaceLabelName(""),
    _buriedEdgesLabelName(""),
    _surfaceLabelValue(1),
    _buriedEdgesLabelValue(1) {
    setLabelValue(100);

    _refDir1[0] = 0.0;
    _refDir1[1] = 0.0;
    _refDir1[2] = 1.0;

    _refDir2[0] = 0.0;
    _refDir2[1] = 1.0;
    _refDir2[2] = 0.0;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesive::~FaultCohesive(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesive::deallocate(void) {
    PYLITH_METHOD_BEGIN;
    pylith::problems::Physics::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
    delete _diagnosticFactory;_diagnosticFactory = NULL;
    delete _derivedFactory;_derivedFactory = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of label identifying cohesive cells.
void
pylith::faults::FaultCohesive::setCohesiveLabelName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setCohesiveLabelName(value="<<value<<")");

    setLabelName(value);
}


// ------------------------------------------------------------------------------------------------
// Get name of label identifying cohesive cells.
const char*
pylith::faults::FaultCohesive::getCohesiveLabelName(void) const {
    return getLabelName();
}


// ------------------------------------------------------------------------------------------------
// Set value of label identifying cohesive cells.
void
pylith::faults::FaultCohesive::setCohesiveLabelValue(const int value) {
    setLabelValue(value);
}


// ------------------------------------------------------------------------------------------------
// Get value of label identifying cohesive cells.
int
pylith::faults::FaultCohesive::getCohesiveLabelValue(void) const {
    return getLabelValue();
}


// ------------------------------------------------------------------------------------------------
// Set name of label marking surface of interface.
void
pylith::faults::FaultCohesive::setSurfaceLabelName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setSurfaceLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for name of label for fault surface.");
    } // if

    _surfaceLabelName = value;
}


// ------------------------------------------------------------------------------------------------
// Get name of label marking surface of interface.
const char*
pylith::faults::FaultCohesive::getSurfaceLabelName(void) const {
    return _surfaceLabelName.c_str();
}


// ------------------------------------------------------------------------------------------------
// Set value of label marking surface of interface.
void
pylith::faults::FaultCohesive::setSurfaceLabelValue(const int value) {
    _surfaceLabelValue = value;
}


// ------------------------------------------------------------------------------------------------
// Get value of label marking surface of interface.
int
pylith::faults::FaultCohesive::getSurfaceLabelValue(void) const {
    return _surfaceLabelValue;
}


// ------------------------------------------------------------------------------------------------
// Set name of label marking buried edges of interface surface.
void
pylith::faults::FaultCohesive::setBuriedEdgesLabelName(const char* value) {
    PYLITH_COMPONENT_DEBUG("setBuriedEdgesLabelName(value="<<value<<")");

    _buriedEdgesLabelName = value;
}


// ------------------------------------------------------------------------------------------------
// Get name of label marking buried edges of interface surface.
const char*
pylith::faults::FaultCohesive::getBuriedEdgesLabelName(void) const {
    return _buriedEdgesLabelName.c_str();
}


// ------------------------------------------------------------------------------------------------
// Set value of label marking buried edges of interface surface.
void
pylith::faults::FaultCohesive::setBuriedEdgesLabelValue(const int value) {
    _buriedEdgesLabelValue = value;
}


// ------------------------------------------------------------------------------------------------
// Get value of label marking buried edges of interface surface.
int
pylith::faults::FaultCohesive::getBuriedEdgesLabelValue(void) const {
    return _buriedEdgesLabelValue;
}


// ------------------------------------------------------------------------------------------------
// Set first choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::faults::FaultCohesive::setRefDir1(const double vec[3]) {
    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 1 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]<<") is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    } // for
} // setRefDir1


// ------------------------------------------------------------------------------------------------
// Set second choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::faults::FaultCohesive::setRefDir2(const double vec[3]) {
    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 2 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]<<") is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir2[i] = vec[i] / mag;
    } // for
} // setRefDir2


// ------------------------------------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(topology::Mesh* const mesh) {
    PYLITH_METHOD_BEGIN;

    assert(mesh);
    assert(_surfaceLabelName.length() > 0);

    try {
        pylith::topology::Mesh faultMesh;

        // Get group of vertices associated with fault
        PetscDM dmMesh = mesh->getDM();assert(dmMesh);

        PetscDMLabel surfaceLabel = NULL;
        PetscBool hasLabel = PETSC_FALSE;
        PetscInt depth, gdepth, dim;
        PetscMPIInt rank;
        PetscErrorCode err;
        // We do not have labels on all ranks until after distribution
        err = MPI_Comm_rank(PetscObjectComm((PetscObject) dmMesh), &rank);PYLITH_CHECK_ERROR(err);
        err = DMHasLabel(dmMesh, _surfaceLabelName.c_str(), &hasLabel);PYLITH_CHECK_ERROR(err);
        if (!hasLabel && !rank) {
            std::ostringstream msg;
            msg << "Mesh missing group of vertices '" << _surfaceLabelName
                << "' for fault interface condition.";
            throw std::runtime_error(msg.str());
        } // if
        err = DMGetDimension(dmMesh, &dim);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetDepth(dmMesh, &depth);PYLITH_CHECK_ERROR(err);
        err = MPI_Allreduce(&depth, &gdepth, 1, MPIU_INT, MPI_MAX, mesh->getComm());PYLITH_CHECK_ERROR(err);
        err = DMGetLabel(dmMesh, _surfaceLabelName.c_str(), &surfaceLabel);PYLITH_CHECK_ERROR(err);
        TopologyOps::createFault(&faultMesh, *mesh, surfaceLabel, _surfaceLabelValue);
        PetscDMLabel buriedEdgesLabel = NULL;

        // We do not have labels on all ranks until after distribution
        if ((_buriedEdgesLabelName.length() > 0) && !rank) {
            err = DMGetLabel(dmMesh, _buriedEdgesLabelName.c_str(), &buriedEdgesLabel);PYLITH_CHECK_ERROR(err);
            if (!buriedEdgesLabel) {
                std::ostringstream msg;
                msg << "Could not find label '" << _buriedEdgesLabelName << "' marking buried edges for fault '" << _surfaceLabelName << "'.";
                throw std::runtime_error(msg.str());
            } // if
        } // if
        TopologyOps::create(mesh, faultMesh, buriedEdgesLabel, _buriedEdgesLabelValue, getCohesiveLabelValue());

        // Check consistency of mesh.
        pylith::topology::MeshOps::checkTopology(*mesh);
        pylith::topology::MeshOps::checkTopology(faultMesh);

        pythia::journal::debug_t debug(PyreComponent::getName());
        if (debug.state()) {
            mesh->view("::ascii_info_detail");
        } // if

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while adjusting topology to create cohesive cells for fault '" << _surfaceLabelName << "'.\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // adjustTopology


// ------------------------------------------------------------------------------------------------
// Create diagnostic field.
pylith::topology::Field*
pylith::faults::FaultCohesive::createDiagnosticField(const pylith::topology::Field& solution,
                                                     const pylith::topology::Mesh& physicsMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDiagnosticField(solution="<<solution.getLabel()<<", physicsMesh=)"<<typeid(physicsMesh).name()<<")");

    assert(_normalizer);

    pylith::topology::Field* diagnosticField = new pylith::topology::Field(physicsMesh);assert(diagnosticField);
    diagnosticField->setLabel("diagnostic field");
    pylith::topology::FieldOps::createOutputLabel(diagnosticField);

    assert(_diagnosticFactory);
    const pylith::topology::FieldBase::Discretization& discretization = solution.getSubfieldInfo("lagrange_multiplier_fault").fe;
    const PylithInt cellDim = solution.getSpaceDim()-1;
    const bool isFaultOnly = false;
    _diagnosticFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder, cellDim,
                                                  isFaultOnly, discretization.cellBasis, discretization.feSpace,
                                                  discretization.isBasisContinuous);

    assert(_diagnosticFactory);
    assert(_normalizer);
    _diagnosticFactory->initialize(diagnosticField, *_normalizer, solution.getSpaceDim());

    _diagnosticFactory->addNormalDir(); // 0
    _diagnosticFactory->addStrikeDir(); // 1
    if (solution.getSpaceDim() > 2) {
        _diagnosticFactory->addUpDipDir(); // 2
    } // if

    diagnosticField->subfieldsSetup();
    diagnosticField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *diagnosticField);
    diagnosticField->allocate();
    diagnosticField->createOutputVector();

    PYLITH_METHOD_RETURN(diagnosticField);
} // createDiagnosticField


// ------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::faults::FaultCohesive::createDerivedField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    assert(_normalizer);

    pylith::topology::Field* derivedField = new pylith::topology::Field(domainMesh);assert(derivedField);
    derivedField->setLabel("derived field");

    // Create label for output.
    const char* outputLabelName = "output";
    PetscDM derivedDM = derivedField->getDM();
    PetscDMLabel outputLabel = NULL;
    PetscErrorCode err = PETSC_SUCCESS;
    err = DMCreateLabel(derivedDM, outputLabelName);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(derivedDM, outputLabelName, &outputLabel);PYLITH_CHECK_ERROR(err);
    pylith::topology::Stratum faultStratum(derivedDM, pylith::topology::Stratum::HEIGHT, 1);
    for (PetscInt point = faultStratum.begin(); point != faultStratum.end(); ++point) {
        err = DMLabelSetValue(outputLabel, point, 1);
    } // for
    err = DMPlexLabelComplete(derivedDM, outputLabel);PYLITH_CHECK_ERROR(err);

    assert(_derivedFactory);
    const pylith::topology::FieldBase::Discretization& discretization = solution.getSubfieldInfo("lagrange_multiplier_fault").fe;
    const PylithInt cellDim = solution.getSpaceDim()-1;
    const bool isFaultOnly = false;
    _derivedFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder, cellDim,
                                               isFaultOnly, discretization.cellBasis, discretization.feSpace,
                                               discretization.isBasisContinuous);

    assert(_derivedFactory);
    assert(_normalizer);
    _derivedFactory->initialize(derivedField, *_normalizer, solution.getSpaceDim());

    _derivedFactory->addTractionChange(); // 0

    derivedField->subfieldsSetup();
    derivedField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *derivedField);
    derivedField->allocate();
    derivedField->createOutputVector();

    PYLITH_METHOD_RETURN(derivedField);
} // createDerivedField


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::faults::FaultCohesive::createIntegrator(const pylith::topology::Field& solution,
                                                const std::vector<pylith::materials::Material*>& materials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorInterface* integrator = new pylith::feassemble::IntegratorInterface(this);assert(integrator);
    integrator->setLabelName(getCohesiveLabelName());
    integrator->setLabelValue(getCohesiveLabelValue());
    integrator->setSurfaceLabelName(getSurfaceLabelName());

    pylith::feassemble::InterfacePatches* patches =
        pylith::feassemble::InterfacePatches::createMaterialPairs(this, solution.getDM());
    integrator->setIntegrationPatches(patches);

    _setKernelsResidual(integrator, solution, materials);
    _setKernelsJacobian(integrator, solution, materials);
    _setKernelsDiagnosticField(integrator, solution);
    _setKernelsDerivedField(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create constraint for buried fault edges and faces.
std::vector<pylith::feassemble::Constraint*>
pylith::faults::FaultCohesive::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getLabel()<<")");

    if (0 == strlen(getBuriedEdgesLabelName())) {
        std::vector<pylith::feassemble::Constraint*> constraintArray;
        PYLITH_METHOD_RETURN(constraintArray);
    } // if

    const char* lagrangeName = "lagrange_multiplier_fault";
    const PylithInt numComponents = solution.getSpaceDim();

    pylith::int_array constrainedDOF;
    constrainedDOF.resize(numComponents);
    for (int c = 0; c < numComponents; ++c) {
        constrainedDOF[c] = c;
    }
    // Make new label for cohesive edges and faces
    PetscDM dm = solution.getDM();
    PetscDMLabel buriedLabel = NULL;
    PetscDMLabel buriedCohesiveLabel = NULL;
    PetscIS pointIS = NULL;
    const PetscInt *points = NULL;
    PetscInt numPoints = 0;
    std::ostringstream labelstream;
    labelstream << getBuriedEdgesLabelName() << "_cohesive";
    std::string buriedLabelName = labelstream.str();
    PetscErrorCode err;

    err = DMCreateLabel(dm, buriedLabelName.c_str());PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, getBuriedEdgesLabelName(), &buriedLabel);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dm, buriedLabelName.c_str(), &buriedCohesiveLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(buriedLabel, getBuriedEdgesLabelValue(), &pointIS);PYLITH_CHECK_ERROR(err);
    if (pointIS) {
        err = ISGetLocalSize(pointIS, &numPoints);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    } // if
    for (int p = 0; p < numPoints; ++p) {
        const PetscInt *support = NULL;
        PetscInt supportSize;

        err = DMPlexGetSupportSize(dm, points[p], &supportSize);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetSupport(dm, points[p], &support);PYLITH_CHECK_ERROR(err);
        for (int s = 0; s < supportSize; ++s) {
            DMPolytopeType ct;
            const PetscInt spoint = support[s];

            err = DMPlexGetCellType(dm, spoint, &ct);PYLITH_CHECK_ERROR(err);
            if ((ct == DM_POLYTOPE_SEG_PRISM_TENSOR) || (ct == DM_POLYTOPE_POINT_PRISM_TENSOR)) {
                const PetscInt *cone = NULL;
                PetscInt coneSize;

                err = DMPlexGetConeSize(dm, spoint, &coneSize);PYLITH_CHECK_ERROR(err);
                err = DMPlexGetCone(dm, spoint, &cone);PYLITH_CHECK_ERROR(err);
                for (int c = 0; c < coneSize; ++c) {
                    PetscInt val;
                    err = DMLabelGetValue(buriedLabel, cone[c], &val);PYLITH_CHECK_ERROR(err);
                    if (val >= 0) {
                        err = DMLabelSetValue(buriedCohesiveLabel, spoint, 1);PYLITH_CHECK_ERROR(err);
                        break;
                    } // if
                } // for
            } // if
        } // for
    } // for
    if (pointIS) {
        err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
    } // if

    std::vector<pylith::feassemble::Constraint*> constraintArray;
    pylith::feassemble::ConstraintSimple *constraint = new pylith::feassemble::ConstraintSimple(this);assert(constraint);
    constraint->setLabelName(buriedLabelName.c_str());
    err = PetscObjectViewFromOptions((PetscObject) buriedLabel, NULL, "-buried_edge_label_view");
    err = PetscObjectViewFromOptions((PetscObject) buriedCohesiveLabel, NULL, "-buried_cohesive_edge_label_view");
    constraint->setConstrainedDOF(&constrainedDOF[0], constrainedDOF.size());
    constraint->setSubfieldName(lagrangeName);
    constraint->setUserFn(_zero);

    constraintArray.resize(1);
    constraintArray[0] = constraint;
    PYLITH_METHOD_RETURN(constraintArray);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::faults::FaultCohesive::_updateKernelConstants(const PylithReal dt) {
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
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::FaultCohesive::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Get derived factory associated with fault.
pylith::topology::FieldFactory*
pylith::faults::FaultCohesive::_getDerivedFactory(void) {
    return _derivedFactory;
} // _getDerivedFactory


// ------------------------------------------------------------------------------------------------
pylith::feassemble::Integrator*
pylith::faults::FaultCohesive::createIntegrator(const pylith::topology::Field& solution) {
    return NULL;
} // Empty method


// ------------------------------------------------------------------------------------------------
// Set kernels for computing diagnostic field.
void
pylith::faults::FaultCohesive::_setKernelsDiagnosticField(pylith::feassemble::IntegratorInterface* integrator,
                                                          const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDiagnosticField(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    const PylithInt spaceDim = solution.getSpaceDim();
    std::vector<ProjectKernels> kernels(spaceDim);
    kernels[0] = ProjectKernels("normal_dir", pylith::fekernels::BoundaryDirections::normalDir);
    kernels[1] = ProjectKernels("strike_dir", pylith::fekernels::BoundaryDirections::tangentialDirHoriz);
    if (spaceDim > 2) {
        kernels[2] = ProjectKernels("up_dip_dir", pylith::fekernels::BoundaryDirections::tangentialDirVert);
    } // if

    assert(integrator);
    integrator->setKernelsDiagnosticField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDiagnosticField


// ------------------------------------------------------------------------------------------------
// Set kernels for computing derived field.
void
pylith::faults::FaultCohesive::_setKernelsDerivedField(pylith::feassemble::IntegratorInterface* integrator,
                                                       const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDerivedField(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels(1);
    kernels[0] = ProjectKernels("traction_change", pylith::fekernels::FaultCohesiveKin::tractionChange_asVector);

    assert(integrator);
    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
