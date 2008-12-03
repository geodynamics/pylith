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

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/CohesiveTopology.hh" // USES CohesiveTopology

#include <Selection.hh> // USES submesh algorithms

#include <map> // USES std::map

// ----------------------------------------------------------------------
// Initialize mesh.
void
pylith::meshio::TestDataWriterVTKBCMesh::_initialize(void)
{ // _initialize
  CPPUNIT_ASSERT(0 != _data);

  MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(&_meshDomain);
  CPPUNIT_ASSERT(!_meshDomain.isNull());

  if (0 != _data->faultLabel) {
    faults::FaultCohesiveKin fault;
    fault.label(_data->faultLabel);
    fault.id(_data->faultId);
    fault.adjustTopology(_meshDomain);
  } // if

  const char* label = _data->bcLabel;
  _mesh = 
    ALE::Selection<Mesh>::submeshV<SubMesh>(_meshDomain, 
					    _meshDomain->getIntSection(label));
  CPPUNIT_ASSERT(!_mesh.isNull());
  _mesh->setRealSection("coordinates", 
			_meshDomain->getRealSection("coordinates"));
  //_mesh->view("BC mesh");
} // _initialize


// End of file 
