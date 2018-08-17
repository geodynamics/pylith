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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesive.hh" // implementation of object methods

#include "pylith/faults/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/faults/TopologyOps.hh" // USES TopologyOps

#include "TopologyOps.hh" // USES TopologyOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDiscretization
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> \
    // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesive::FaultCohesive(void) :
    _faultMesh(NULL),
    _cohesivePointMap(NULL),
    _auxFaultFactory(new pylith::faults::AuxiliaryFactory),
    _id(100),
    _label(""),
    _edge(""){ // constructor
    _refDir1[0] = 0.0;
    _refDir1[1] = 0.0;
    _refDir1[2] = 1.0;

    _refDir2[0] = 0.0;
    _refDir2[1] = 1.0;
    _refDir2[2] = 0.0;
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesive::~FaultCohesive(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesive::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::feassemble::IntegratorPointwise::deallocate();
    ISDestroy(&_cohesivePointMap);
    delete _faultMesh;_faultMesh = NULL;
    delete _auxFaultFactory;_auxFaultFactory = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set material identifier of fault.
void
pylith::faults::FaultCohesive::id(const int value) {
    PYLITH_COMPONENT_DEBUG("id(value="<<value<<")");

    _id = value;
} // id


// ----------------------------------------------------------------------
// Get material identifier of fault.
int
pylith::faults::FaultCohesive::id(void) const {
    return _id;
} // id


// ----------------------------------------------------------------------
// Set label of group of vertices associated with fault.
void
pylith::faults::FaultCohesive::label(const char* value) {
    PYLITH_COMPONENT_DEBUG("label(value="<<value<<")");

    _label = value;
} // label


// ----------------------------------------------------------------------
// Get label of group of vertices associated with fault.
const char*
pylith::faults::FaultCohesive::label(void) const {
    return _label.c_str();
} // label


// ----------------------------------------------------------------------
// Set label of group of vertices defining buried edge of fault.
void
pylith::faults::FaultCohesive::edge(const char* value) {
    PYLITH_COMPONENT_DEBUG("edge(value="<<value<<")");

    _edge = value;
} // edge


// ----------------------------------------------------------------------
// Get label of group of vertices defining buried edge of fault.
const char*
pylith::faults::FaultCohesive::edge(void) const {
    return _edge.c_str();
} // edge


// ----------------------------------------------------------------------
// Set first choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::faults::FaultCohesive::refDir1(const double vec[3]) {
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
} // refDir1


// ----------------------------------------------------------------------
// Set second choice for reference direction to discriminate among tangential directions in 3-D.
void
pylith::faults::FaultCohesive::refDir2(const double vec[3]) {
    // Set reference direction, insuring it is a unit vector.
    const PylithReal mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (mag < 1.0e-6) {
        std::ostringstream msg;
        msg << "Magnitude of reference direction 2 ("<<vec[0]<<", "<<vec[1]<<", "<<vec[2]<<") is negligible. Use a unit vector.";
        throw std::runtime_error(msg.str());
    } // if
    for (int i = 0; i < 3; ++i) {
        _refDir1[i] = vec[i] / mag;
    } // for
} // refDir2


// ----------------------------------------------------------------------
// Get mesh associated with integrator domain.
const pylith::topology::Mesh&
pylith::faults::FaultCohesive::domainMesh(void) const {
    assert(!_faultMesh);
    return *_faultMesh;
} // domainMesh


// ----------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(topology::Mesh* const mesh){ // adjustTopology
    PYLITH_METHOD_BEGIN;

    assert(mesh);
    assert(std::string("") != label());

    try {
        pylith::topology::Mesh faultMesh;

        // Get group of vertices associated with fault
        PetscDM dmMesh = mesh->dmMesh();assert(dmMesh);

        const char* charlabel = label();

        PetscDMLabel groupField;
        PetscBool hasLabel;
        PetscInt depth, gdepth, dim;
        PetscMPIInt rank;
        PetscErrorCode err;
        // We do not have labels on all ranks until after distribution
        err = MPI_Comm_rank(PetscObjectComm((PetscObject) dmMesh), &rank);PYLITH_CHECK_ERROR(err);
        err = DMHasLabel(dmMesh, charlabel, &hasLabel);PYLITH_CHECK_ERROR(err);
        if (!hasLabel && !rank) {
            std::ostringstream msg;
            msg << "Mesh missing group of vertices '" << label()
                << "' for fault interface condition.";
            throw std::runtime_error(msg.str());
        } // if
        err = DMGetDimension(dmMesh, &dim);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetDepth(dmMesh, &depth);PYLITH_CHECK_ERROR(err);
        err = MPI_Allreduce(&depth, &gdepth, 1, MPIU_INT, MPI_MAX, mesh->comm());PYLITH_CHECK_ERROR(err);
        err = DMGetLabel(dmMesh, charlabel, &groupField);PYLITH_CHECK_ERROR(err);
        TopologyOps::createFault(&faultMesh, *mesh, groupField);
        PetscDMLabel faultBdLabel = NULL;

        // We do not have labels on all ranks until after distribution
        if (( strlen(edge()) > 0) && !rank) {
            err = DMGetLabel(dmMesh, edge(), &faultBdLabel);PYLITH_CHECK_ERROR(err);
            if (!faultBdLabel) {
                std::ostringstream msg;
                msg << "Could not find nodeset/pset '" << edge() << "' marking buried edges for fault '" << label() << "'.";
                throw std::runtime_error(msg.str());
            } // if
        } // if
        TopologyOps::create(mesh, faultMesh, faultBdLabel, id());

        // Check consistency of mesh.
        pylith::topology::MeshOps::checkTopology(*mesh);
        pylith::topology::MeshOps::checkTopology(faultMesh);

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while adjusting topology to create cohesive cells for fault '" << label() << "'.\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // adjustTopology


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesive::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize integrator. Create fault mesh from cohesive cells and cohesive point map.
void
pylith::faults::FaultCohesive::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label()<<")");

    const bool isSubMesh = true;
    delete _faultMesh;_faultMesh = new pylith::topology::Mesh(isSubMesh);assert(_faultMesh);
    pylith::faults::TopologyOps::createFaultParallel(_faultMesh, solution.mesh(), id(), label());
    pylith::topology::MeshOps::checkTopology(*_faultMesh);

    // Optimize closure for coordinates.
    PetscDM dmFault = _faultMesh->dmMesh();assert(dmFault);
    pylith::topology::CoordsVisitor::optimizeClosure(dmFault);

    // Set default discretization of auxiliary subfields to match lagrange_multiplier_fault subfield in solution.
    assert(_auxFaultFactory);
    const pylith::topology::FieldBase::Discretization& discretization = solution.subfieldInfo("lagrange_multiplier_fault").fe;
    _auxFaultFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder,
                                                discretization.isBasisContinuous, discretization.feSpace);

    delete _auxField;_auxField = new pylith::topology::Field(*_faultMesh);assert(_auxField);
    _auxField->label("fault auxiliary");
    _auxFieldSetup();
    _auxField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *_auxField);
    _auxField->allocate();
    _auxField->zeroLocal();

    assert(_normalizer);
    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory();assert(factory);
    factory->initializeSubfields();

    PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Get factory for setting up auxliary fields.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::FaultCohesive::_auxFactory(void) {
    return _auxFaultFactory;
} // _auxFactory


// End of file
