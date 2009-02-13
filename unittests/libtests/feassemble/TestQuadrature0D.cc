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
#include "pylith/feassemble/QuadratureRefCell.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestQuadrature0D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestQuadrature0D::testConstructor(void)
{ // testConstructor
  QuadratureRefCell refCell;
  Quadrature0D quadrature(refCell);
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
  
  QuadratureRefCell refCell;
  refCell.minJacobian(minJacobian);
  refCell.initialize(basis, basisDeriv, quadPtsRef, quadWts,
		     cellDim, numBasis, numQuadPts, spaceDim);

  Quadrature0D engine(refCell);

  engine.initialize();
  engine.computeGeometry(vertCoords, spaceDim, 0);

  const double tolerance = 1.0e-06;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(quadPts[0], engine._quadPts[0], 
			       tolerance);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobian[0], engine._jacobian[0], 
			       tolerance);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianInv[0], engine._jacobianInv[0], 
			       tolerance);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(jacobianDet[0], engine._jacobianDet[0], 
			       tolerance);
} // testPoint


// End of file 
