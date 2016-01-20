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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultMesh.hh" // Implementation of class methods

#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology
#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
void
pylith::faults::TestFaultMesh::createFaultMesh(topology::Mesh* faultMesh,
					       topology::Mesh* mesh,
					       const char* faultLabel,
					       const int faultId)
{ // createFaultMesh
  PYLITH_METHOD_BEGIN;
  
  CPPUNIT_ASSERT(faultMesh);
  CPPUNIT_ASSERT(mesh);

  PetscErrorCode err = 0;
  { // Create mesh
    PetscInt depth = -1, firstFaultVertex = 0;
    PetscInt firstLagrangeVertex = 0, firstFaultCell = 0;
    PetscDMLabel groupField = NULL;
    const bool useLagrangeConstraints = true;
    PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
    
    err = DMGetStratumSize(dmMesh, faultLabel, 1, &firstLagrangeVertex);PYLITH_CHECK_ERROR(err);
    firstFaultCell = firstLagrangeVertex;
    if (useLagrangeConstraints) {
      firstFaultCell += firstLagrangeVertex;
    } // if
    err = DMPlexGetDepth(dmMesh, &depth);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmMesh, faultLabel, &groupField);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT(groupField);
    CohesiveTopology::createFault(faultMesh, *mesh, groupField);
    CohesiveTopology::create(mesh, *faultMesh, NULL, faultId, firstFaultVertex, firstLagrangeVertex, firstFaultCell, useLagrangeConstraints);
  } // Create mesh

  PYLITH_METHOD_END;
} // createFaultMesh
