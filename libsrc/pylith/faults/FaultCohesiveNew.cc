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

#include "FaultCohesiveNew.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveNew::FaultCohesiveNew(void) :
    _faultMesh(NULL),
    _cohesivePointMap(NULL),
    _id(100)
{ // constructor
    _label = "";
    _edge = "";
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveNew::~FaultCohesiveNew(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesiveNew::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::feassemble::IntegratorPointwise::deallocate();
    delete _faultMesh; _faultMesh = NULL;
    ISDestroy(_cohesivePointMap);

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set material identifier of fault.
void
pylith::faults::FaultCohesiveNew::id(const int value) {
    PYLITH_COMPONENT_DEBUG("id(value="<<value<<")");

    _id = value;
} // id

// ----------------------------------------------------------------------
// Get material identifier of fault.
int
pylith::faults::FaultCohesiveNew::id(void) const {
    return _id;
} // id


// ----------------------------------------------------------------------
// Set label of group of vertices associated with fault.
void
pylith::faults::FaultCohesiveNew::label(const char* value) {
    PYLITH_COMPONENT_DEBUG("label(value="<<value<<")");

    _label = value;
} // label

// ----------------------------------------------------------------------
// Get label of group of vertices associated with fault.
const char*
pylith::faults::FaultCohesiveNew::label(void) const {
    return _label.c_str();
} // label

// ----------------------------------------------------------------------
// Set label of group of vertices defining buried edge of fault.
void
pylith::faults::FaultCohesiveNew::edge(const char* value) {
    PYLITH_COMPONENT_DEBUG("edge(value="<<value<<")");

    _edge = value;
} // edge

// ----------------------------------------------------------------------
// Get label of group of vertices defining buried edge of fault.
const char*
pylith::faults::FaultCohesiveNew::edge(void) const {
    return _edge.c_str();
} // edge

// ----------------------------------------------------------------------
// Set up direction to discriminate among shear directions in 3-D.
void
pylith::faults::FaultCohesiveNew::upDir(const double vec[3]) {
    PYLITH_COMPONENT_DEBUG("upDir(vec=["<<vec[0]<<","<<vec[1]<<","<<vec[2]<<"])");

    for (int i = 0; i < 3; ++i) {
        _upDir[i] = vec[i];
    } // for
} // upDir


// ----------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesiveNew::adjustTopology(topology::Mesh* const mesh)
{ // adjustTopology
    PYLITH_METHOD_BEGIN;

    assert(mesh);
    assert(std::string("") != label());

    try {
        topology::Mesh faultMesh;

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
        CohesiveTopology::createFault(&faultMesh, *mesh, groupField);
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
        CohesiveTopology::create(mesh, faultMesh, faultBdLabel, id(), *firstFaultVertex, *firstLagrangeVertex, *firstFaultCell, useLagrangeConstraints());

        // Check consistency of mesh.
        topology::MeshOps::checkTopology(*mesh);
        topology::MeshOps::checkTopology(faultMesh);

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while adjusting topology to create cohesive cells for fault '" << label() << "'.\n"
            << err.what();
        throw std::runtime_error(msg.str());
    }

    PYLITH_METHOD_END;
} // adjustTopology


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesiveNew::verifyConfiguration(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label() < ")");

    PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Initialize integrator. Create fault mesh from cohesive cells and cohesive point map.
void
pylith::faults::FaultCohesiveNew::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(solution="<<solution.label() < ")");

    PYLITH_METHOD_END;
} // initialize


// End of file
