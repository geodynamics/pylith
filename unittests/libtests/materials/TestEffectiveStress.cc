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

#include "TestEffectiveStress.hh" // Implementation of class methods

#include "pylith/materials/EffectiveStress.hh" // USES EffectiveStress

#include <cmath> // USES pow()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestEffectiveStress );

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _EffectiveStress {
      class Linear {
      public :
	Linear(void) {};
	~Linear(void) {};
	double effStressFunc(const double x) {
	  return x - 10.0;
	};
	double effStressDerivFunc(const double x) {
	  return 1.0;
	};
	double effStressFuncDerivFunc(double* f,
				      double* df,
				      const double x) {
	  *f = effStressFunc(x);
	  *df = effStressDerivFunc(x);
	};
      }; // Linear
    } // _EffectiveStress
  } // materials
} // pylith

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _EffectiveStress {
      class Quadratic {
      public :
	Quadratic(void) {};
	~Quadratic(void) {};
	double effStressFunc(const double x) {
	  return 1.0e+5 - 1.0/9.0e+3 * pow(x + 2.0e+4, 2);
	};
	double effStressDerivFunc(const double x) {
	  return -2*1.0/9.0e+3*(x+2.0e+4);
	};
	double effStressFuncDerivFunc(double* f,
				      double* df,
				      const double x) {
	  *f = effStressFunc(x);
	  *df = effStressDerivFunc(x);
	};
      }; // Quadratic
    } // _EffectiveStress
  } // materials
} // pylith

// ----------------------------------------------------------------------
namespace pylith {
  namespace materials {
    namespace _EffectiveStress {
      class Cubic {
      public :
	Cubic(void) {};
	~Cubic(void) {};
	double effStressFunc(const double x) {
	  return pow(x - 4.0, 3) - 8.0;
	};
	double effStressDerivFunc(const double x) {
	  return 3.0*pow(x - 4.0, 2);
	};
	double effStressFuncDerivFunc(double* f,
				      double* df,
				      const double x) {
	  *f = effStressFunc(x);
	  *df = effStressDerivFunc(x);
	};
      }; // Cubic
    } // _EffectiveStress
  } // materials
} // pylith

// ----------------------------------------------------------------------
// Test calculate() with linear function.
void
pylith::materials::TestEffectiveStress::testCalculateLinear(void)
{ // testCalculateLinear
  const double valueE = 10.0;
  
  _EffectiveStress::Linear material;

  const int ntests = 4;
  const double guesses[ntests] = { 0.0, 6.0, 14.0, 20.0 };
  const double scale = 1.0;
  const double tolerance = 1.0e-06;
  for (int i=0; i < ntests; ++i) {
    const double value =
      EffectiveStress::calculate<_EffectiveStress::Linear>(guesses[i], scale,
							   &material);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, value/valueE, tolerance);
  } // for
} // testCalculateLinear

// ----------------------------------------------------------------------
// Test calculate() with quadratic function.
void
pylith::materials::TestEffectiveStress::testCalculateQuadratic(void)
{ // testCalculateQuadratic
  const double valueE = 1.0e+04;
  
  _EffectiveStress::Quadratic material;

  const int ntests = 4;
  const double guesses[ntests] = { 1.0, 1.0e-1, 2.0e-2, 1.0e-2 };
  const double scale = 1.0e-2;
  const double tolerance = 1.0e-06;
  for (int i=0; i < ntests; ++i) {
    const double value =
      EffectiveStress::calculate<_EffectiveStress::Quadratic>(guesses[i], scale,
							   &material);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, value/valueE, tolerance);
  } // for
} // testCalculateQuadratic

// ----------------------------------------------------------------------
// Test calculate() with cubic function.
void
pylith::materials::TestEffectiveStress::testCalculateCubic(void)
{ // testCalculateCubic
  const double valueE = 6.0;
  
  _EffectiveStress::Cubic material;

  const int ntests = 4;
  const double guesses[ntests] = { 2.0, 4.0, 6.0, 8.0 };
  const double scale = 1.0;
  const double tolerance = 1.0e-06;
  for (int i=0; i < ntests; ++i) {
    const double value =
      EffectiveStress::calculate<_EffectiveStress::Cubic>(guesses[i], scale,
							   &material);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, value/valueE, tolerance);
  } // for
} // testCalculateCubic


// End of file
