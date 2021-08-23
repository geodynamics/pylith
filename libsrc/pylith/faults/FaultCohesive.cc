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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesive.hh" // implementation of object methods

#include "pylith/faults/TopologyOps.hh" // USES TopologyOps
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/FEKernelKey.hh" // USES FEKernelKey
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/MeshOps.hh" // USES MeshOps::checkTopology()

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <utility> // USES std::pair
#include <map> // USES std::map
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesive::FaultCohesive(void) :
    _interfaceId(100),
    _interfaceLabel(""),
    _buriedEdgesLabel("") {
    _refDir1[0] = 0.0;
    _refDir1[1] = 0.0;
    _refDir1[2] = 1.0;

    _refDir2[0] = 0.0;
    _refDir2[1] = 1.0;
    _refDir2[2] = 0.0;
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesive::~FaultCohesive(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesive::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::problems::Physics::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set identifier for fault cohesive cells.
void
pylith::faults::FaultCohesive::setInterfaceId(const int value) {
    PYLITH_COMPONENT_DEBUG("setInterfaceId(value="<<value<<")");

    _interfaceId = value;
} // setInterfaceId


// ---------------------------------------------------------------------------------------------------------------------
// Get identifier for fault cohesive cells.
int
pylith::faults::FaultCohesive::getInterfaceId(void) const {
    return _interfaceId;
} // getInterfaceId


// ---------------------------------------------------------------------------------------------------------------------
// Set label marking surface of interface.
void
pylith::faults::FaultCohesive::setSurfaceMarkerLabel(const char* value) {
    PYLITH_COMPONENT_DEBUG("setSurfaceMarkerLabel(value="<<value<<")");

    _interfaceLabel = value;
} // setSurfaceMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label marking surface of interface.
const char*
pylith::faults::FaultCohesive::getSurfaceMarkerLabel(void) const {
    return _interfaceLabel.c_str();
} // getSurfaceMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Set label marking buried edges of interface surface.
void
pylith::faults::FaultCohesive::setBuriedEdgesMarkerLabel(const char* value) {
    PYLITH_COMPONENT_DEBUG("setBuriedEdgesMarkerLabel(value="<<value<<")");

    _buriedEdgesLabel = value;
} // setBuriedEdgesMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label marking buried edges of interface surface.
const char*
pylith::faults::FaultCohesive::getBuriedEdgesMarkerLabel(void) const {
    return _buriedEdgesLabel.c_str();
} // getBuriedEdgesMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
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


// ---------------------------------------------------------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(topology::Mesh* const mesh) {
    PYLITH_METHOD_BEGIN;

    assert(mesh);
    assert(_interfaceLabel.length() > 0);

    try {
        pylith::topology::Mesh faultMesh;

        // Get group of vertices associated with fault
        PetscDM dmMesh = mesh->getDM();assert(dmMesh);

        PetscDMLabel groupField = NULL;
        PetscBool hasLabel = PETSC_FALSE;
        PetscInt depth, gdepth, dim;
        PetscMPIInt rank;
        PetscErrorCode err;
        // We do not have labels on all ranks until after distribution
        err = MPI_Comm_rank(PetscObjectComm((PetscObject) dmMesh), &rank);PYLITH_CHECK_ERROR(err);
        err = DMHasLabel(dmMesh, _interfaceLabel.c_str(), &hasLabel);PYLITH_CHECK_ERROR(err);
        if (!hasLabel && !rank) {
            std::ostringstream msg;
            msg << "Mesh missing group of vertices '" << _interfaceLabel
                << "' for fault interface condition.";
            throw std::runtime_error(msg.str());
        } // if
        err = DMGetDimension(dmMesh, &dim);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetDepth(dmMesh, &depth);PYLITH_CHECK_ERROR(err);
        err = MPI_Allreduce(&depth, &gdepth, 1, MPIU_INT, MPI_MAX, mesh->getComm());PYLITH_CHECK_ERROR(err);
        err = DMGetLabel(dmMesh, _interfaceLabel.c_str(), &groupField);PYLITH_CHECK_ERROR(err);
        TopologyOps::createFault(&faultMesh, *mesh, groupField);
        PetscDMLabel faultBdLabel = NULL;

        // We do not have labels on all ranks until after distribution
        if ((_buriedEdgesLabel.length() > 0) && !rank) {
            err = DMGetLabel(dmMesh, _buriedEdgesLabel.c_str(), &faultBdLabel);PYLITH_CHECK_ERROR(err);
            if (!faultBdLabel) {
                std::ostringstream msg;
                msg << "Could not find nodeset/pset '" << _buriedEdgesLabel << "' marking buried edges for fault '" << _interfaceLabel << "'.";
                throw std::runtime_error(msg.str());
            } // if
        } // if
        TopologyOps::create(mesh, faultMesh, faultBdLabel, _interfaceId);

        // Check consistency of mesh.
        pylith::topology::MeshOps::checkTopology(*mesh);
        pylith::topology::MeshOps::checkTopology(faultMesh);

    } catch (const std::exception& err) {
        std::ostringstream msg;
        msg << "Error occurred while adjusting topology to create cohesive cells for fault '" << _interfaceLabel << "'.\n"
            << err.what();
        throw std::runtime_error(msg.str());
    } // try/catch

    PYLITH_METHOD_END;
} // adjustTopology


// End of file
