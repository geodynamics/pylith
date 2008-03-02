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

#include "TestDataWriterVTKMesh.hh" // Implementation of class methods

#include "data/DataWriterVTKData.hh" // USES DataWriterVTKData

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterVTKMesh::_initialize(void)
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);

  MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(&_mesh);
  CPPUNIT_ASSERT(!_mesh.isNull());
  _mesh->getFactory()->clear();

  if (0 != _data->faultLabel) {
    faults::FaultCohesiveKin fault;
    fault.label(_data->faultLabel);
    fault.id(_data->faultId);
    fault.adjustTopology(_mesh);
  } // if
} // _initialize


// End of file 
