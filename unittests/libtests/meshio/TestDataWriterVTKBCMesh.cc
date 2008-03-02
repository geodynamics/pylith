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

#include "TestDataWriterVTKBCMesh.hh" // Implementation of class methods

#include "data/DataWriterVTKData.hh" // USES DataWriterVTKData
#include "data/DataWriterVTKDataMesh.hh" // USES DataWriterVTKDataMesh

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology

#include <Selection.hh> // USES submesh algorithms

#include <map> // USES std::map

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKBCMesh::setUp(void)
{ // setUp
  TestDataWriterVTK::setUp();
  _dataMesh = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKBCMesh::tearDown(void)
{ // tearDown
  TestDataWriterVTK::tearDown();
  delete _dataMesh; _dataMesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterVTKBCMesh::_initialize(void)
{ // _initialize
  CPPUNIT_ASSERT(0 != _dataMesh);

  MeshIOAscii iohandler;
  iohandler.filename(_dataMesh->meshFilename);
  iohandler.read(&_meshDomain);
  CPPUNIT_ASSERT(!_meshDomain.isNull());
  _meshDomain->getFactory()->clear();

  if (0 != _dataMesh->faultLabel) {
    faults::FaultCohesiveKin fault;
    fault.label(_dataMesh->faultLabel);
    fault.id(_dataMesh->faultId);
    fault.adjustTopology(_meshDomain);
  } // if

  const char* label = _dataMesh->bcLabel;
  _mesh = 
    ALE::Selection<ALE::Mesh>::submesh(_meshDomain, 
				       _meshDomain->getIntSection(label));
  CPPUNIT_ASSERT(!_mesh.isNull());
  _mesh->setRealSection("coordinates", 
			_meshDomain->getRealSection("coordinates"));
  //_mesh->view("BC mesh");
} // _initialize


// End of file 
