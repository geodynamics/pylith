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

#include "TestFaultCohesive.hh" // Implementation of class methods

//#include "data/FaultCohesiveData.hh" // USES FaultCohesiveData

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultsCohesiveKin

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

#include <stdexcept> // TEMPORARY

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesive );

// ----------------------------------------------------------------------
// Test adjustTopology() with 1-D line element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyLine(void)
{ // testAdjustTopologyLine
  const char* filename = "data/meshTet4A_orig.txt";
  _testAdjustTopology(filename);
} // testAdjustTopologyLine

// ----------------------------------------------------------------------
// Test adjustTopology().
void
pylith::faults::TestFaultCohesive::_testAdjustTopology(const char* filename)
{ // _testAdjustTopology
  ALE::Obj<ALE::Mesh> mesh;
  meshio::MeshIOAscii iohandler;
  iohandler.filename(filename);
  iohandler.debug(true);
  iohandler.interpolate(true);
  iohandler.read(&mesh);

  FaultCohesiveKin fault;
  fault.id(0);
  fault.label("fault");
  fault.adjustTopology(&mesh);

  mesh->view("Mesh");

  throw std::logic_error("Unit test not fully implemented.");
} // _testAdjustTopology

// End of file 
