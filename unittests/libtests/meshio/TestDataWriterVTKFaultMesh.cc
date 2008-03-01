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

#include "TestDataWriterVTKFaultMesh.hh" // Implementation of class methods

#include "data/DataWriterVTKData.hh" // USES DataWriterVTKData
#include "data/DataWriterVTKDataMesh.hh" // USES DataWriterVTKDataMesh

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology

#include <map> // USES std::map

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTKFaultMesh::setUp(void)
{ // setUp
  TestDataWriterVTK::setUp();
  _dataMesh = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTKFaultMesh::tearDown(void)
{ // tearDown
  TestDataWriterVTK::tearDown();
  delete _dataMesh; _dataMesh = 0;
} // tearDown

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterVTKFaultMesh::_initialize(void)
{ // _initialize
  CPPUNIT_ASSERT(0 != _dataMesh);

  MeshIOAscii iohandler;
  iohandler.filename(_dataMesh->meshFilename);
  iohandler.read(&_meshDomain);
  CPPUNIT_ASSERT(!_meshDomain.isNull());
  _meshDomain->getFactory()->clear();

  faults::FaultCohesiveKin fault;
  fault.label(_dataMesh->faultLabel);
  fault.id(_dataMesh->faultId);
  fault.adjustTopology(_meshDomain);
  const bool constraintCell = true;
  std::map<Mesh::point_type, Mesh::point_type> cohesiveToFault;
  faults::CohesiveTopology::createParallel(&_mesh, &cohesiveToFault,
					   _meshDomain, _dataMesh->faultId,
					   constraintCell);
} // _initialize


// End of file 
