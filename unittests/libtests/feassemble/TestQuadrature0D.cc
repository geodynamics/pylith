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

#include "TestQuadrature0D.hh" // Implementation of class methods

#include "pylith/feassemble/Quadrature0D.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature0D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature0D::testConstructor(void)
{ // testConstructor
  Quadrature0D quadrature;
} // testConstructor

// ----------------------------------------------------------------------
// Test computeGeometry().
void
pylith::feassemble::TestQuadrature0D::testPoint(void)
{ // testPoint
  const int cellDim = 0;
  const int numBasis = 1;
  const int numQuadPts = 1;
  const int spaceDim = 1;
  const double verticesRef[] = { 0.0 };
  const double basis[] = { 1.0 };
  const double basisDeriv[] = { 1.0 };
  const double quadPtsRef[] = { 0.0 };
  const double quadWts[] = { 1.0 };

  const int numVertices = 1;
  const int numCells = 1;
  const double vertCoords[] = { 1.1 };
  const int cells[] = { 0 };
  const double quadPts[] = { 1.1 };
  const double jacobian[] = { 1.0 };
  const double jacobianInv[] = { 1.0 };
  const double jacobianDet[] = { 1.0 };

  const double minJacobian = 1.0e-06;
  
  Quadrature0D quadrature;

  quadrature.minJacobian(minJacobian);
  quadrature.initialize(verticesRef, basis, basisDeriv, quadPtsRef, quadWts,
		    cellDim, numBasis, numQuadPts, spaceDim);

  // Create mesh with test cell
  ALE::Obj<Mesh> mesh = new Mesh(PETSC_COMM_WORLD, cellDim);
  CPPUNIT_ASSERT(!mesh.isNull());
  ALE::Obj<sieve_type> sieve = new sieve_type(mesh->comm());
  CPPUNIT_ASSERT(!sieve.isNull());

  const bool interpolate = false;
  ALE::SieveBuilder<Mesh>::buildTopology(sieve, cellDim, numCells,
		     (int*) cells, numVertices, interpolate, numBasis);
  mesh->setSieve(sieve);
  mesh->stratify();
  ALE::SieveBuilder<Mesh>::buildCoordinates(mesh, spaceDim, vertCoords);
  
  // Check values from computeGeometry()
  const ALE::Obj<Mesh::label_sequence>& cellsMesh = mesh->heightStratum(0);
  CPPUNIT_ASSERT(!cellsMesh.isNull());
  const Mesh::label_sequence::iterator e_iter = cellsMesh->begin(); 
  const ALE::Obj<Mesh::real_section_type>& coordinates = 
    mesh->getRealSection("coordinates");
  CPPUNIT_ASSERT(!coordinates.isNull());
  quadrature.computeGeometry(mesh, coordinates, *e_iter);

  CPPUNIT_ASSERT(1 == numCells);

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPts[0], quadrature._quadPts[0], 
			       tolerance);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobian[0], quadrature._jacobian[0], 
			       tolerance);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianInv[0], quadrature._jacobianInv[0], 
			       tolerance);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDet[0], quadrature._jacobianDet[0], 
			       tolerance);
} // testPoint


// End of file 
