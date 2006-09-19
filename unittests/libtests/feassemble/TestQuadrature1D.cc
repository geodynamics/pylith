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

#include "TestQuadrature1D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature1D.hh" // USES Quadrature1D

#include <sstream> // USES std::stringstream

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature1D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature1D::testConstructor(void)
{ // testConstructor
  Quadrature1D quadrature;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry() w/linear basis fns
void
pylith::feassemble::TestQuadrature1D::testLinear(void)
{ // testLinear
  // Basis fns: N0 = 1-p, N1 = p
  // Quadrature points: (p=0.5, wt=1.0)

  const int cellDim = 1;
  const int numCorners = 2;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double basis[] = { 0.5, 0.5 };
  const double basisDeriv[] = { -1.0, 1.0 };
  const double quadPtsRef[] = { 0.5 };
  const double quadWts[] = { 1.0 };
  const double jacobianTol = 1.0;

  Quadrature1D q;
  q.initialize(basis, basisDeriv, quadPtsRef, quadWts,
	       cellDim, numCorners, numQuadPts, spaceDim);

  // Test cell: x0=-0.25, x1=2.0 
  const int numVertices = 2;
  const int numCells = 1;
  const double vertCoords[] = { -0.25, 2.0 };
  const int cells[] = { 0, 1 };
  const double quadPts[] = { 0.875 };
  const double jacobian[] = { 2.25 };
  const double jacobianInv[] = { 1.0/2.25 };
  const double jacobianDet[] = { 2.25 };

  // Create mesh with test cell
  typedef ALE::Sieve<int, int, int> sieve_type;
  typedef ALE::New::Topology<int, sieve_type> topology_type;
  typedef ALE::New::Section<topology_type, double> section_type;
  ALE::Obj<ALE::Mesh> mesh = ALE::Mesh(PETSC_COMM_WORLD, cellDim);
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());
  ALE::Obj<topology_type> topology = new topology_type(mesh->comm());

  const bool interpolate = false;
  ALE::New::SieveBuilder<sieve_type>::buildTopology(sieve, cellDim, numCells,
		     (int*) cells, numVertices, interpolate, numCorners);
  sieve->stratify();
  topology->setPatch(0, sieve);
  topology->stratify();
  mesh->setTopologyNew(topology);
  ALE::New::SieveBuilder<sieve_type>::buildCoordinates(
		    mesh->getSection("coordinates"), spaceDim, vertCoords);
  
  // Check values from _computeGeometry()
  const ALE::Mesh::topology_type::patch_type patch = 0;
  const ALE::Obj<topology_type::label_sequence>& elements = 
    topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator e_iter = elements->begin(); 
  const ALE::Obj<ALE::Mesh::section_type>& coordinates = 
    mesh->getSection("coordinates");
  q._computeGeometry(coordinates, *e_iter);

  int size = 1;
  CPPUNIT_ASSERT(0 != q._quadPts);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(quadPts[i], q._quadPts[i]);

  size = 1;
  CPPUNIT_ASSERT(0 != q._jacobian);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobian[i], q._jacobian[i]);

  size = 1;
  CPPUNIT_ASSERT(0 != q._jacobianInv);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianInv[i], q._jacobianInv[i]);

  size = 1;
  CPPUNIT_ASSERT(0 != q._jacobianDet);
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_EQUAL(jacobianDet[i], q._jacobianDet[i]);
} // testLinear

// ----------------------------------------------------------------------
// Test computeGeometry() w/quadratic basis fns
void
pylith::feassemble::TestQuadrature1D::testQuadratic(void)
{ // testQuadratic
  CPPUNIT_ASSERT(false);
} // testQuadratic

// End of file 
