// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIntegratorElasticity.hh" // Implementation of class methods

#include "pylith/feassemble/IntegratorElasticity.hh" // USES IntegratorElasticity

#include <math.h> // USES fabs()

#include <stdexcept>
// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestIntegratorElasticity );

// ----------------------------------------------------------------------
// Test calcTotalStrain1D().
void
pylith::feassemble::TestIntegratorElasticity::testCalcTotalStrain1D(void)
{ // testCalcTotalStrain1D
  // N0 = 0.5 * (1 - x)
  // N1 = 0.5 * (1 + x)
  // dN0/dx = -0.5
  // dN1/dx = +0.5
  // Let quad pt 0 be dN/dx, let quad pt 1 be 0.5*dN/dx
  const int dim = 1;
  const int numBasis = 2;
  const int numQuadPts = 2;
  const PylithScalar basisDerivVals[] = {
    -0.50, 0.50,
    -0.25, 0.25 };
  const int tensorSize = 1;

  // Let u(x) = 1 + 0.5 * x
  const PylithScalar dispVals[] = { 0.5, 1.5 };
  const PylithScalar strainE[] = { 0.5, 0.25 };

  const int size = numQuadPts * tensorSize;
  scalar_array strain(size);

  scalar_array basisDeriv(basisDerivVals, numQuadPts*numBasis*dim);
  scalar_array disp(dispVals, numBasis*dim);

  IntegratorElasticity::_calcTotalStrain1D(&strain,
					   basisDeriv, disp, numBasis, numQuadPts);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_EQUAL(size, int(strain.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(strainE[i], strain[i], tolerance);
} // testCalcTotalStrain1D

// ----------------------------------------------------------------------
// Test calcTotalStrain2D().
void
pylith::feassemble::TestIntegratorElasticity::testCalcTotalStrain2D(void)
{ // testCalcTotalStrain2D
  // N0 = x
  // N1 = y
  // N2 = 1 - x - y
  // dN0/dx = +1.0, dN0/dy =  0.0
  // dN1/dx =  0.0, dN1/dy = +1.0
  // dN2/dx = -1.0, dN2/dy = -1.0
  // Let quad pt 0 be derivatives, let quad pt 1 be 2.0*derivatives
  const int dim = 2;
  const int numBasis = 3;
  const int numQuadPts = 2;
  const PylithScalar basisDerivVals[] = {
    +1.0,  0.0,   0.0, +1.0,   -1.0, -1.0,
    +2.0,  0.0,   0.0, +2.0,   -2.0, -2.0
  };
  const int tensorSize = 3;

  // Let ux(x,y) = +0.4 + 0.3*x + 0.8*y
  // Ley uy(x,y) = -2.0 + 0.5*x - 0.2*y
  const PylithScalar dispVals[] = {
    0.7, -1.5,
    1.2, -2.2,
    0.4, -2.0
  };
  const PylithScalar strainE[] = {
    0.3, -0.2, 0.65,
    0.6, -0.4, 1.3
  };

  const int size = numQuadPts * tensorSize;
  scalar_array strain(size);

  scalar_array basisDeriv(basisDerivVals, numQuadPts*numBasis*dim);
  scalar_array disp(dispVals, numBasis*dim);

  IntegratorElasticity::_calcTotalStrain2D(&strain, basisDeriv, disp,
					   numBasis, numQuadPts);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_EQUAL(size, int(strain.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(strainE[i], strain[i], tolerance);
} // testCalcTotalStrain2D

// ----------------------------------------------------------------------
// Test calcTotalStrain3D().
void
pylith::feassemble::TestIntegratorElasticity::testCalcTotalStrain3D(void)
{ // testCalcTotalStrain3D
  // N0 = x
  // N1 = y
  // N2 = z
  // N2 = 1 - x - y - z
  // dN0/dx = +1.0, dN0/dy =  0.0, dN0/dz =  0.0
  // dN1/dx =  0.0, dN1/dy = +1.0, dN0/dz =  0.0
  // dN2/dx =  0.0, dN2/dy =  0.0, dN2/dz = +1.0
  // dN3/dx = -1.0, dN3/dy = -1.0, dN3/dz = -1.0
  // Let quad pt 0 be derivatives, let quad pt 1 be 3.0*derivatives
  const int dim = 3;
  const int numBasis = 4;
  const int numQuadPts = 2;
  const PylithScalar basisDerivVals[] = {
    +1.0,  0.0,  0.0, // Quad pt 0
     0.0, +1.0,  0.0,
     0.0,  0.0, +1.0,
    -1.0, -1.0, -1.0,
    +3.0,  0.0,  0.0, // Quad pt 1
     0.0, +3.0,  0.0,
     0.0,  0.0, +3.0,
    -3.0, -3.0, -3.0
  };
  const int tensorSize = 6;

  // Let ux(x,y,z) = +0.4 + 0.3*x + 0.8*y + 0.4*z
  // Ley uy(x,y,z) = -2.0 + 0.5*x - 0.2*y + 1.2*z
  // Ley uz(x,y,z) = -1.0 + 0.2*x - 0.7*y - 0.3*z
  const PylithScalar dispVals[] = {
    0.7, -1.5, -0.8,
    1.2, -2.2, -1.7,
    0.8, -0.8, -1.3,
    0.4, -2.0, -1.0
  };
  const PylithScalar strainE[] = {
    0.3, -0.2, -0.3, 0.65, 0.25, 0.3,
    0.9, -0.6, -0.9, 1.95, 0.75, 0.9
  };

  const int size = numQuadPts * tensorSize;
  scalar_array strain(size);

  scalar_array basisDeriv(basisDerivVals, numQuadPts*numBasis*dim);
  scalar_array disp(dispVals, numBasis*dim);

  IntegratorElasticity::_calcTotalStrain3D(&strain, basisDeriv, disp,
					   numBasis, numQuadPts);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_EQUAL(size, int(strain.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(strainE[i], strain[i], tolerance);
} // testCalcTotalStrain3D


// End of file 
