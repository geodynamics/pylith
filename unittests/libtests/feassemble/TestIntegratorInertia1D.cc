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

#include "TestIntegratorInertia1D.hh" // Implementation of class methods

#include "pylith/feassemble/IntegratorInertia.hh"
#include "pylith/feassemble/Quadrature1D.hh"

#include "data/IntegratorDataInertia1DLinear.hh"
#include "data/IntegratorDataInertia1DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestIntegratorInertia1D );

// ----------------------------------------------------------------------
// Test constructor
void
pylith::feassemble::TestIntegratorInertia1D::testConstructor(void)
{ // testConstructor
  IntegratorInertia integrator;
} // testConstructor

// ----------------------------------------------------------------------
// Test integrateAction() w/linear basis fns
void
pylith::feassemble::TestIntegratorInertia1D::testActionLinear(void)
{ // testActionLinear
  IntegratorDataInertia1DLinear data;

  Quadrature1D q;
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
// Test integrate() w/linear basis fns
void
pylith::feassemble::TestIntegratorInertia1D::testIntegrateLinear(void)
{ // testIntegrateLinear
  IntegratorDataInertia1DLinear data;

  Quadrature1D q;
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
// Test integrateLumped() w/linear basis fns
void
pylith::feassemble::TestIntegratorInertia1D::testLumpedLinear(void)
{ // testLumpedLinear
  IntegratorDataInertia1DLinear data;

  Quadrature1D q;
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
// Test integrateAction() w/quadratic basis fns
void
pylith::feassemble::TestIntegratorInertia1D::testActionQuadratic(void)
{ // testActionQuadratic
  IntegratorDataInertia1DQuadratic data;

  Quadrature1D q;
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
} // testActionQuadratic

// ----------------------------------------------------------------------
// Test integrate() w/quadratic basis fns
void
pylith::feassemble::TestIntegratorInertia1D::testIntegrateQuadratic(void)
{ // testIntegrateQuadratic
  IntegratorDataInertia1DQuadratic data;

  Quadrature1D q;
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
} // testIntegrateQuadratic

// ----------------------------------------------------------------------
// Test integrateLumped() w/quadratic basis fns
void
pylith::feassemble::TestIntegratorInertia1D::testLumpedQuadratic(void)
{ // testLumpedQuadratic
  IntegratorDataInertia1DQuadratic data;

  Quadrature1D q;
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
} // testLumpedQuadratic

// End of file 
