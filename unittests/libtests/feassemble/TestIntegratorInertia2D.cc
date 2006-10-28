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

#include "TestIntegratorInertia2D.hh" // Implementation of class methods

#include "pylith/feassemble/IntegratorInertia.hh"
#include "pylith/feassemble/Quadrature2Din3D.hh"

#include "data/IntegratorDataInertia2Din3DOne.hh"
#include "data/IntegratorDataInertia2Din3DTwo.hh"
#include "data/IntegratorDataInertia2Din3DThree.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestIntegratorInertia2D );

// ----------------------------------------------------------------------
// Test integrate() & integrateAction() w/linear basis fns (1 cell)
void
pylith::feassemble::TestIntegratorInertia2D::testOne(void)
{ // testOne
  IntegratorDataInertia2Din3DOne data;

  Quadrature2Din3D q;
  q.initialize(data.basis,
	       data.basisDeriv,
	       data.quadPts,
	       data.quadWts,
	       data.cellDim,
	       data.numCorners,
	       data.numQuadPts,
	       data.spaceDim);

  IntegratorInertia integrator;
  integrator.quadrature(&q);

  _testIntegrateAction(&integrator, data);

  _testIntegrate(&integrator, data);
} // testLinear

// ----------------------------------------------------------------------
// Test integrate() & integrateAction() w/linear basis fns
// (2 cells sharing 1 vertex)
void
pylith::feassemble::TestIntegratorInertia2D::testOverlap1(void)
{ // testOverlap1
  IntegratorDataInertia2Din3DTwo data;
  
  Quadrature2Din3D q;
  q.initialize(data.basis,
	       data.basisDeriv,
	       data.quadPts,
	       data.quadWts,
	       data.cellDim,
	       data.numCorners,
	       data.numQuadPts,
	       data.spaceDim);

  IntegratorInertia integrator;
  integrator.quadrature(&q);

  _testIntegrateAction(&integrator, data);

  _testIntegrate(&integrator, data);
} // testOverlap1

// ----------------------------------------------------------------------
// Test integrate() & integrateAction() w/linear basis fns
// (2 cells sharing 2 vertices)
void
pylith::feassemble::TestIntegratorInertia2D::testOverlap2(void)
{ // testOverlap2
  IntegratorDataInertia2Din3DThree data;
  
  Quadrature2Din3D q;
  q.initialize(data.basis,
	       data.basisDeriv,
	       data.quadPts,
	       data.quadWts,
	       data.cellDim,
	       data.numCorners,
	       data.numQuadPts,
	       data.spaceDim);

  IntegratorInertia integrator;
  integrator.quadrature(&q);

  _testIntegrateAction(&integrator, data);

  _testIntegrate(&integrator, data);
} // testOverlap2

// End of file 
