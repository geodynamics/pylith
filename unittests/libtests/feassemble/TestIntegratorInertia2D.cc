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
// Test integrateAction() w/linear basis fns (1 cell)
void
pylith::feassemble::TestIntegratorInertia2D::testActionOne(void)
{ // testActionOne
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
} // testActionLinear

// ----------------------------------------------------------------------
// Test integrate() w/linear basis fns (1 cell)
void
pylith::feassemble::TestIntegratorInertia2D::testIntegrateOne(void)
{ // testIntegrateOne
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

  _testIntegrate(&integrator, data);
} // testIntegrateLinear

// ----------------------------------------------------------------------
// Test integrateLumped() w/linear basis fns (1 cell)
void
pylith::feassemble::TestIntegratorInertia2D::testLumpedOne(void)
{ // testLumpedOne
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

  _testIntegrateLumped(&integrator, data);
} // testLumpedLinear

// ----------------------------------------------------------------------
// Test integrateAction() w/linear basis fns (2 cells sharing 1 vertex)
void
pylith::feassemble::TestIntegratorInertia2D::testActionOverlap1(void)
{ // testActionOverlap1
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
} // testActionOverlap1

// ----------------------------------------------------------------------
// Test integrateAction() w/linear basis fns (2 cells sharing 1 vertex)
void
pylith::feassemble::TestIntegratorInertia2D::testIntegrateOverlap1(void)
{ // testIntegrateOverlap1
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

  _testIntegrate(&integrator, data);
} // testIntegrateOverlap1

// ----------------------------------------------------------------------
// Test integrateLumped() w/linear basis fns (2 cells sharing 1 vertex)
void
pylith::feassemble::TestIntegratorInertia2D::testLumpedOverlap1(void)
{ // testLumpedOverlap1
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

  _testIntegrateLumped(&integrator, data);
} // testLumpedOverlap1

// ----------------------------------------------------------------------
// Test integrateAction() w/linear basis fns (2 cells sharing 2 vertices)
void
pylith::feassemble::TestIntegratorInertia2D::testActionOverlap2(void)
{ // testActionOverlap2
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
} // testActionOverlap2

// ----------------------------------------------------------------------
// Test integrate() w/linear basis fns (2 cells sharing 2 vertices)
void
pylith::feassemble::TestIntegratorInertia2D::testIntegrateOverlap2(void)
{ // testIntegrateOverlap2
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

  _testIntegrate(&integrator, data);
} // testIntegrateOverlap2

// ----------------------------------------------------------------------
// Test integrateLumped() w/linear basis fns (2 cells sharing 2 vertices)
void
pylith::feassemble::TestIntegratorInertia2D::testLumpedOverlap2(void)
{ // testLumpedOverlap2
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

  _testIntegrateLumped(&integrator, data);
} // testLumpedOverlap2

// End of file 
