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

#include "TestSolutionFields.hh" // Implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::topology::TestSolutionFields );

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestSolutionFields::testConstructor(void)
{ // testConstructor
  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);
} // testConstructor
 
// ----------------------------------------------------------------------
// Test solutionName().
void
pylith::topology::TestSolutionFields::testSolutionName(void)
{ // testSolutionName
  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);

  const std::string& name = "my solution";
  manager.add(name.c_str(), "displacement");
  manager.solutionName(name.c_str());
  CPPUNIT_ASSERT_EQUAL(name, manager._solutionName);
} // testSolutionName

// ----------------------------------------------------------------------
// Test solution().
void
pylith::topology::TestSolutionFields::testSolution(void)
{ // testSolution
  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);

  const char* labels[] = { "field A", "field B", "field C" };
  const int size = 3;
  const int fiberDimA = 2;
  const int fiberDimB = 3;
  const int fiberDimC = 4;

  for (int i=0; i < size; ++i)
    manager.add(labels[i], "displacement");

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  Field<Mesh>& fieldA = manager.get(labels[0]);
  Field<Mesh>& fieldB = manager.get(labels[1]);
  Field<Mesh>& fieldC = manager.get(labels[2]);
  fieldA.newSection(vertices, fiberDimA);

  const ALE::Obj<Mesh::RealSection>& section = fieldA.section();
  const Mesh::RealSection::chart_type& chart = section->getChart();
  fieldB.newSection(chart, fiberDimB);
  fieldC.newSection(chart, fiberDimC);

  manager.solutionName(labels[1]);
  const Field<Mesh>& solution = manager.solution();
  const ALE::Obj<Mesh::RealSection>& sectionSoln = solution.section();
  CPPUNIT_ASSERT_EQUAL(fiberDimB,
		       sectionSoln->getFiberDimension(*(vertices->begin())));
} // testSolution

// ----------------------------------------------------------------------
// Test createHistory().
void
pylith::topology::TestSolutionFields::testCreateHistory(void)
{ // testCreateHistory
  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);

  const char* labels[] = { "field A", "field B", "field C" };
  const int totalSize = 3;
  const int historySize = 2;

  // Add fields
  for (int i=0; i < totalSize; ++i)
    manager.add(labels[i], "displacement");

  manager.createHistory(labels, historySize);
  for (int i=0; i < historySize; ++i)
    CPPUNIT_ASSERT_EQUAL(std::string(labels[i]), manager._history[i]);
} // testCreateHistory

// ----------------------------------------------------------------------
// Test shiftHistory().
void
pylith::topology::TestSolutionFields::testShiftHistory(void)
{ // testShiftHistory
  Mesh mesh;
  _initialize(&mesh);
  SolutionFields manager(mesh);

  const char* fieldNames[] = { "field A", "field B" };
  const int numFields = 2;
  const int fiberDimA = 2;
  const int fiberDimB = 3;

  for (int i=0; i < numFields; ++i)
    manager.add(fieldNames[i], "displacement");
  manager.createHistory(fieldNames, numFields);

  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = mesh.sieveMesh();
  const ALE::Obj<Mesh::SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  Field<Mesh>& fieldA = manager.get(fieldNames[0]);
  Field<Mesh>& fieldB = manager.get(fieldNames[1]);
  fieldA.newSection(vertices, fiberDimA);
  fieldB.newSection(vertices, fiberDimB);

  manager.shiftHistory();
  const Field<Mesh>& testA = manager.get(fieldNames[0]);
  const ALE::Obj<Mesh::RealSection>& sectionA = testA.section();
  const Field<Mesh>& testB = manager.get(fieldNames[1]);
  const ALE::Obj<Mesh::RealSection>& sectionB = testB.section();
  CPPUNIT_ASSERT_EQUAL(fiberDimB, 
		       sectionA->getFiberDimension(*(vertices->begin())));
  CPPUNIT_ASSERT_EQUAL(fiberDimA, 
		       sectionB->getFiberDimension(*(vertices->begin())));

} // testShiftHistory

// ----------------------------------------------------------------------
void
pylith::topology::TestSolutionFields::_initialize(Mesh* mesh) const
{ // _initialize
  CPPUNIT_ASSERT(0 != mesh);

  meshio::MeshIOAscii iohandler;
  iohandler.filename("data/tri3.mesh");
  iohandler.read(mesh);
} // _initialize


// End of file 
