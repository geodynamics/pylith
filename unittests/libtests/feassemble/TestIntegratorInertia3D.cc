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

#include "TestIntegratorInertia3D.hh" // Implementation of class methods

#include "pylith/feassemble/IntegratorInertia.hh"
#include "pylith/feassemble/Quadrature3D.hh"

#include "data/IntegratorDataInertia3DLinear.hh"
#include "data/IntegratorDataInertia3DQuadratic.hh"

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestIntegratorInertia3D );

// ----------------------------------------------------------------------
// Test integrateAction() w/linear basis fns
void
pylith::feassemble::TestIntegratorInertia3D::testActionLinear(void)
{ // testActionLinear
  IntegratorDataInertia3DLinear data;

  Quadrature3D q;
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
pylith::feassemble::TestIntegratorInertia3D::testIntegrateLinear(void)
{ // testIntegrateLinear
  IntegratorDataInertia3DLinear data;

  Quadrature3D q;
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
pylith::feassemble::TestIntegratorInertia3D::testLumpedLinear(void)
{ // testLumpedLinear
  IntegratorDataInertia3DLinear data;

  Quadrature3D q;
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
pylith::feassemble::TestIntegratorInertia3D::testActionQuadratic(void)
{ // testActionQuadratic
  IntegratorDataInertia3DQuadratic data;

  Quadrature3D q;
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
pylith::feassemble::TestIntegratorInertia3D::testIntegrateQuadratic(void)
{ // testIntegrateQuadratic
  IntegratorDataInertia3DQuadratic data;

  Quadrature3D q;
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
pylith::feassemble::TestIntegratorInertia3D::testLumpedQuadratic(void)
{ // testLumpedQuadratic
  IntegratorDataInertia3DQuadratic data;

  Quadrature3D q;
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
