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
// Test integrate() & integrateAction() w/linear basis fns
void
pylith::feassemble::TestIntegratorInertia3D::testLinear(void)
{ // testLinear
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

  _testIntegrate(&integrator, data);
} // testLinear

// ----------------------------------------------------------------------
// Test integrate() & integrateAction() w/quadratic basis fns
void
pylith::feassemble::TestIntegratorInertia3D::testQuadratic(void)
{ // testQuadratic
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

  _testIntegrate(&integrator, data);
} // testQuadratic

// End of file 
