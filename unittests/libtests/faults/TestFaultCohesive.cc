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

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesive );

// ----------------------------------------------------------------------
// Test adjustTopology() with 1-D line element.
void
pylith::faults::TestFaultCohesive::testAdjustTopologyLine(void)
{ // testAdjustTopologyLine
  //FaultDataLine data;
  //_testAdjustTopologyLine(data);
} // testAdjustTopologyLine

// ----------------------------------------------------------------------
// Test adjustTopology().
void
pylith::faults::TestFaultCohesive::_testAdjustTopologyLine(
						      const FaultData& data)
{ // _testAdjustTopologyLine
  //ALE::Obj<ALE::Mesh> mesh(new ALE::Mesh);
  //_createMesh(&mesh, data);

  FaultCohesiveKin fault;
} // _testAdjustTopologyLine

// ----------------------------------------------------------------------
// Create mesh.
void
pylith::faults::TestFaultCohesive::_createMesh(ALE::Obj<ALE::Mesh>* mesh,
					       const FaultData& data)
{ // _createMesh
} // _createMesh


// End of file 
