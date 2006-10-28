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
// Test integrate() & integrateAction() w/linear basis fns
void
pylith::feassemble::TestIntegratorInertia1D::testLinear(void)
{ // testLinear
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

  _testIntegrate(&integrator, data);
} // testLinear

// ----------------------------------------------------------------------
// Test integrate() & integrateAction() w/quadratic basis fns
void
pylith::feassemble::TestIntegratorInertia1D::testQuadratic(void)
{ // testQuadratic
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

  _testIntegrate(&integrator, data);
} // testQuadratic

// End of file 
