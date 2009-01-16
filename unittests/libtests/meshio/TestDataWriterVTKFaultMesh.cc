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

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology

#include <map> // USES std::map

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterVTKFaultMesh::_initialize(void)
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);

  MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(&_meshDomain);
  CPPUNIT_ASSERT(!_meshDomain.isNull());

  faults::FaultCohesiveKin fault;
  fault.label(_data->faultLabel);
  fault.id(_data->faultId);
  fault.adjustTopology(_meshDomain, _flipFault);
  const bool constraintCell = true;
  std::map<Mesh::point_type, Mesh::point_type> cohesiveToFault;
  faults::CohesiveTopology::createParallel(&_mesh, &cohesiveToFault,
					   _meshDomain, _data->faultId,
					   constraintCell);
} // _initialize


// End of file 
