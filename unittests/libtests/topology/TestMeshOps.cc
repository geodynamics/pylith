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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMeshOps.hh" // Implementation of class methods

#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestMeshOps );

// ----------------------------------------------------------------------
// Test checkMaterialIds().
void
pylith::topology::TestMeshOps::testCheckMaterialIds(void)
{ // testCheckMaterialIds
  Mesh mesh;

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(&mesh);

  const int numMaterials = 2;
  int materialIds[numMaterials];
  materialIds[0] = 4;
  materialIds[1] = 3;

  MeshOps::checkMaterialIds(mesh, materialIds, numMaterials);

  materialIds[0] = 99;
    
  CPPUNIT_ASSERT_THROW(MeshOps::checkMaterialIds(mesh, materialIds, numMaterials),
		       std::runtime_error);
} // testCheckMaterialIds
 

// End of file 
