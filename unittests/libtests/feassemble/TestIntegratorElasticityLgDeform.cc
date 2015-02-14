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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIntegratorElasticityLgDeform.hh" // Implementation of class methods

#include "pylith/feassemble/IntegratorElasticityLgDeform.hh" // USES IntegratorElasticityLgDeform

#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

#include <math.h> // USES fabs()

#include <stdexcept>

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::feassemble::TestIntegratorElasticityLgDeform );

// ----------------------------------------------------------------------
// Test calcDeformation() for 2-D.
void
pylith::feassemble::TestIntegratorElasticityLgDeform::testCalcDeformation2D(void)
{ // testCalcDeformation2D
  PYLITH_METHOD_BEGIN;

  // N0 = x
  // N1 = y
  // N2 = 1 - x - y
  // dN0/dx = +1.0, dN0/dy =  0.0
  // dN1/dx =  0.0, dN1/dy = +1.0
  // dN2/dx = -1.0, dN2/dy = -1.0
  const int dim = 2;
  const int numBasis = 3;
  const int numQuadPts = 2;
  const PylithScalar vertices[numBasis*dim] = {
    1.0, 0.0,
    0.0, 1.0, 
    0.0, 0.0,
  };
  const PylithScalar basisDerivVals[numQuadPts*numBasis*dim] = {
    +1.0,  0.0,   0.0, +1.0,   -1.0, -1.0,
    +1.0,  0.0,   0.0, +1.0,   -1.0, -1.0
  };
  const int tensorSize = 3;

  const int size = numQuadPts * dim*dim;
  scalar_array deform(size);    
  scalar_array basisDeriv(basisDerivVals, numQuadPts*numBasis*dim);

  const PylithScalar tolerance = 1.0e-06;

  { // Rigid body translation
    // ux(x,y) = 0.5
    // uy(x,y) = 0.2
    const PylithScalar disp[numBasis*dim] = {
      0.5, 0.2,
      0.5, 0.2,
      0.5, 0.2,
    };
    const PylithScalar deformE[size] = {
      1.0, 0.0,   0.0, 1.0,
      1.0, 0.0,   0.0, 1.0,
    };
    
    IntegratorElasticityLgDeform::_calcDeformation(&deform, basisDeriv, vertices, disp, numBasis, numQuadPts, dim);
    
    CPPUNIT_ASSERT_EQUAL(size, int(deform.size()));
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(deformE[i], deform[i], tolerance);
  } // Rigid body translation

  { // Rigid body translation + rotation
    // ux(x,y) = 0.5 + cos(theta)*x + sin(theta)*y - x0
    // uy(x,y) = 0.2 - sin(theta)*x + cos(theta)*y - y0
    // theta = pi/6
    const PylithScalar pi = 4.0*atan(1.0);
    const PylithScalar theta = pi / 6.0;
    const PylithScalar disp[numBasis*dim] = {
      0.5+cos(theta)*1.0+sin(theta)*0.0-1.0,
      0.2-sin(theta)*1.0+cos(theta)*0.0-0.0,

      0.5+cos(theta)*0.0+sin(theta)*1.0-0.0,
      0.2-sin(theta)*0.0+cos(theta)*1.0-1.0,

      0.5+cos(theta)*0.0+sin(theta)*0.0-0.0,
      0.2-sin(theta)*0.0+cos(theta)*0.0-0.0,
    };
    const PylithScalar deformE[size] = {
      cos(theta), sin(theta), -sin(theta), cos(theta),
      cos(theta), sin(theta), -sin(theta), cos(theta),
    };
    
    IntegratorElasticityLgDeform::_calcDeformation(&deform, basisDeriv, vertices, disp, numBasis, numQuadPts, dim);
    
    CPPUNIT_ASSERT_EQUAL(size, int(deform.size()));
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(deformE[i], deform[i], tolerance);
  } // Rigid body translation + rotation

  { // Uniform strain
    // Let ux(x,y) = +0.4 + 0.3*x + 0.8*y
    // Ley uy(x,y) = -2.0 + 0.5*x - 0.2*y
    const PylithScalar disp[numBasis*dim] = {
      0.7, -1.5,
      1.2, -2.2,
      0.4, -2.0
    };
    const PylithScalar deformE[size] = {
      1.3, 0.8,   0.5, 0.8,
      1.3, 0.8,   0.5, 0.8,
    };
    
    IntegratorElasticityLgDeform::_calcDeformation(&deform, basisDeriv, vertices, disp, numBasis, numQuadPts, dim);
    
    CPPUNIT_ASSERT_EQUAL(size, int(deform.size()));
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(deformE[i], deform[i], tolerance);
  } // Uniform strain

  PYLITH_METHOD_END;
} // testCalcDeformation2D

// ----------------------------------------------------------------------
// Test calcDeformation() for 3-D.
void
pylith::feassemble::TestIntegratorElasticityLgDeform::testCalcDeformation3D(void)
{ // testCalcDeformation3D
  PYLITH_METHOD_BEGIN;

  // N0 = 0.125 * (1-x) * (1-y) * (1-z)
  // N1 = 0.125 * (1+x) * (1-y) * (1-z)
  // N2 = 0.125 * (1+x) * (1+y) * (1-z)
  // N3 = 0.125 * (1-x) * (1+y) * (1-z)
  // N4 = 0.125 * (1-x) * (1-y) * (1+z)
  // N5 = 0.125 * (1+x) * (1-y) * (1+z)
  // N6 = 0.125 * (1+x) * (1+y) * (1+z)
  // N7 = 0.125 * (1-x) * (1+y) * (1+z)
  // dN0/dx=-0.125, dN0/dy=-0.125, dN0/dz=-0.125
  // dN1/dx=+0.125, dN1/dy=-0.125, dN1/dz=-0.125
  // dN2/dx=+0.125, dN2/dy=+0.125, dN2/dz=-0.125
  // dN3/dx=-0.125, dN3/dy=+0.125, dN3/dz=-0.125
  // dN4/dx=-0.125, dN0/dy=-0.125, dN0/dz=+0.125
  // dN5/dx=+0.125, dN1/dy=-0.125, dN1/dz=+0.125
  // dN6/dx=+0.125, dN2/dy=+0.125, dN2/dz=+0.125
  // dN7/dx=-0.125, dN3/dy=+0.125, dN3/dz=+0.125
  const int dim = 3;
  const int numBasis = 8;
  const int numQuadPts = 1;
  const PylithScalar vertices[numBasis*dim] = {
    -1.0, -1.0, -1.0,
    +1.0, -1.0, -1.0,
    +1.0, +1.0, -1.0,
    -1.0, +1.0, -1.0,
    -1.0, -1.0, +1.0,
    +1.0, -1.0, +1.0,
    +1.0, +1.0, +1.0,
    -1.0, +1.0, +1.0,
  };
  const PylithScalar basisDerivVals[numQuadPts*numBasis*dim] = {
    -0.125, -0.125, -0.125,
    +0.125, -0.125, -0.125,
    +0.125, +0.125, -0.125,
    -0.125, +0.125, -0.125,
    -0.125, -0.125, +0.125,
    +0.125, -0.125, +0.125,
    +0.125, +0.125, +0.125,
    -0.125, +0.125, +0.125,
  };
  const int tensorSize = 3;

  const int size = numQuadPts * dim*dim;
  scalar_array deform(size);    
  scalar_array basisDeriv(basisDerivVals, numQuadPts*numBasis*dim);

  const PylithScalar tolerance = 1.0e-06;

  { // Rigid body translation
    // ux(x,y,z) = 0.5
    // uy(x,y,z) = 0.2
    // uz(x,y,z) = 0.3
    const PylithScalar disp[numBasis*dim] = {
      0.5, 0.2, 0.3,
      0.5, 0.2, 0.3,
      0.5, 0.2, 0.3,
      0.5, 0.2, 0.3,
      0.5, 0.2, 0.3,
      0.5, 0.2, 0.3,
      0.5, 0.2, 0.3,
      0.5, 0.2, 0.3,
    };
    const PylithScalar deformE[size] = {
      1.0, 0.0, 0.0,   0.0, 1.0, 0.0,  0.0, 0.0, 1.0,
    };
    
    IntegratorElasticityLgDeform::_calcDeformation(&deform, basisDeriv, vertices, disp, numBasis, numQuadPts, dim);
    
    CPPUNIT_ASSERT_EQUAL(size, int(deform.size()));
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(deformE[i], deform[i], tolerance);
  } // Rigid body translation

  { // Uniform strain
    // Let ux(x,y,z) = +0.4 + 0.3*x + 0.8*y - 0.2*z
    // Ley uy(x,y,z) = -2.0 + 0.5*x - 0.2*y + 0.6*z
    // Ley uz(x,y,z) = +0.7 + 0.8*x - 0.9*y - 0.1*z
    const PylithScalar disp[numBasis*dim] = {
      +0.4+0.3*-1+0.8*-1-0.2*-1,
      -2.0+0.5*-1-0.2*-1+0.6*-1,
      +0.7+0.8*-1-0.9*-1-0.1*-1,

      +0.4+0.3*+1+0.8*-1-0.2*-1,
      -2.0+0.5*+1-0.2*-1+0.6*-1,
      +0.7+0.8*+1-0.9*-1-0.1*-1,

      +0.4+0.3*+1+0.8*+1-0.2*-1,
      -2.0+0.5*+1-0.2*+1+0.6*-1,
      +0.7+0.8*+1-0.9*+1-0.1*-1,

      +0.4+0.3*-1+0.8*+1-0.2*-1,
      -2.0+0.5*-1-0.2*+1+0.6*-1,
      +0.7+0.8*-1-0.9*+1-0.1*-1,

      +0.4+0.3*-1+0.8*-1-0.2*+1,
      -2.0+0.5*-1-0.2*-1+0.6*+1,
      +0.7+0.8*-1-0.9*-1-0.1*+1,

      +0.4+0.3*+1+0.8*-1-0.2*+1,
      -2.0+0.5*+1-0.2*-1+0.6*+1,
      +0.7+0.8*+1-0.9*-1-0.1*+1,

      +0.4+0.3*+1+0.8*+1-0.2*+1,
      -2.0+0.5*+1-0.2*+1+0.6*+1,
      +0.7+0.8*+1-0.9*+1-0.1*+1,

      +0.4+0.3*-1+0.8*+1-0.2*+1,
      -2.0+0.5*-1-0.2*+1+0.6*+1,
      +0.7+0.8*-1-0.9*+1-0.1*+1,

    };
    const PylithScalar deformE[size] = {
      1.3, 0.8, -0.2,
      0.5, 0.8, 0.6,
      0.8, -0.9, 0.9,
    };
    
    IntegratorElasticityLgDeform::_calcDeformation(&deform, basisDeriv, vertices, disp, numBasis, numQuadPts, dim);
    
    CPPUNIT_ASSERT_EQUAL(size, int(deform.size()));
    for (int i=0; i < size; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(deformE[i], deform[i], tolerance);
  } // Uniform strain

  PYLITH_METHOD_END;
} // testCalcDeformation3D

// ----------------------------------------------------------------------
// Test calcTotalStrain2D().
void
pylith::feassemble::TestIntegratorElasticityLgDeform::testCalcTotalStrain2D(void)
{ // testCalcTotalStrain2D
  PYLITH_METHOD_BEGIN;

  // Deformation tensor X
  // X = [ 2.0, 0.4 ]    [  1.5, 0.2 ]
  //     [ 0.6, 1.4 ],   [ -0.9, 0.8 ]
  const int dim = 2;
  const int numQuadPts = 2;
  const PylithScalar deformVals[] = {
    2.0,  0.4, 0.6, 1.4,
    1.5, 0.2, -0.9, 0.8,
  };
  const PylithScalar strainE[] = {
    1.68, 0.56, 0.82,
    1.03, -0.16, -0.21,
  };
  const int tensorSize = 3;

  const int size = numQuadPts * tensorSize;
  scalar_array strain(size);

  scalar_array deform(deformVals, numQuadPts*dim*dim);

  IntegratorElasticityLgDeform::_calcTotalStrain2D(&strain, deform, numQuadPts);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_EQUAL(size, int(strain.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(strainE[i], strain[i], tolerance);

  PYLITH_METHOD_END;
} // testCalcTotalStrain2D

// ----------------------------------------------------------------------
// Test calcTotalStrain3D().
void
pylith::feassemble::TestIntegratorElasticityLgDeform::testCalcTotalStrain3D(void)
{ // testCalcTotalStrain3D
  PYLITH_METHOD_BEGIN;

  // Deformation tensor X
  // X = [  2.0, 0.4,  0.3 ]    [  1.5, 0.2, -0.1 ]
  //     [  0.6, 1.4, -0.8 ]    [ -0.9, 0.8,  0.3 ]
  //     [ -0.1, 0.6,  1.0 ],   [ -0.8, 0.5,  1.2 ]
  const int dim = 3;
  const int numQuadPts = 2;
  const PylithScalar deformVals[] = {
    2.0,  0.4, 0.3,  0.6, 1.4, -0.8,  -0.1, 0.6, 1.0,
    1.5, 0.2, -0.1,  -0.9, 0.8, 0.3,  -0.8, 0.5, 1.2,
  };
  const PylithScalar strainE[] = {
    1.685, 0.74, 0.365, 0.79, -0.2, 0.01,
    1.35, -0.035, 0.27, -0.41, 0.41, -0.69,
  };
  const int tensorSize = 6;

  const int size = numQuadPts * tensorSize;
  scalar_array strain(size);

  scalar_array deform(deformVals, numQuadPts*dim*dim);

  IntegratorElasticityLgDeform::_calcTotalStrain3D(&strain, deform, numQuadPts);

  const PylithScalar tolerance = 1.0e-06;
  CPPUNIT_ASSERT_EQUAL(size, int(strain.size()));
  for (int i=0; i < size; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(strainE[i], strain[i], tolerance);

  PYLITH_METHOD_END;
} // testCalcTotalStrain3D


// End of file 
