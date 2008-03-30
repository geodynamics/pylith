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

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestMeshOps );

// ----------------------------------------------------------------------
// Test checkMaterialIds().
void
pylith::topology::TestMeshOps::testCheckMaterialIds(void)
{ // testCheckMaterialIds
  ALE::Obj<Mesh> mesh;

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(&mesh);
  CPPUNIT_ASSERT(!mesh.isNull());

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
