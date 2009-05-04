// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestMeshOps.hh" // Implementation of class methods

#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

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

  bool caughtError = false;
  try {
    materialIds[0] = 99;
    
    MeshOps::checkMaterialIds(mesh, materialIds, numMaterials);
  } catch (const std::runtime_error& err) {
    caughtError = true;
  } // try/catch
  CPPUNIT_ASSERT(caughtError);
} // testCheckMaterialIds
 

// End of file 
