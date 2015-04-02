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

/**
 * @file unittests/libtests/faults/TestFaultCohesiveKinCases.hh
 *
 * @brief C++ unit testing for FaultCohesiveKin.
 */

#if !defined(pylith_faults_testfaultcohesivekincases_hh)
#define pylith_faults_testfaultcohesivekincases_hh

#include "TestFaultCohesiveKin.hh" // ISA TestFaultCohesiveKin

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinTri3;
    class TestFaultCohesiveKinTri3d;
    class TestFaultCohesiveKinTri3g;
    class TestFaultCohesiveKinQuad4;
    class TestFaultCohesiveKinQuad4e;
    class TestFaultCohesiveKinQuad4i;
    class TestFaultCohesiveKinTet4;
    class TestFaultCohesiveKinTet4e;
    class TestFaultCohesiveKinTet4f;
    class TestFaultCohesiveKinHex8;
  } // bc
} // pylith

// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 2-D triangular cells.
class pylith::faults::TestFaultCohesiveKinTri3 : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTri3

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTri3 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testAdjustSolnLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinTri3


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 2-D triangular cells.
class pylith::faults::TestFaultCohesiveKinTri3d : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTri3d

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTri3d );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinTri3d


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 2-D triangular cells.
class pylith::faults::TestFaultCohesiveKinTri3g : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTri3g

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTri3g );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinTri3g


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveKinQuad4 : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinQuad4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinQuad4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testAdjustSolnLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinQuad4


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveKinQuad4e : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinQuad4e

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinQuad4e );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinQuad4e

// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 2-D quadrilateral cells.
class pylith::faults::TestFaultCohesiveKinQuad4i : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinQuad4i

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinQuad4i );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinQuad4i


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
class pylith::faults::TestFaultCohesiveKinTet4 : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTet4

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTet4 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testAdjustSolnLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinTet4


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
class pylith::faults::TestFaultCohesiveKinTet4e : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTet4e

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTet4e );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinTet4e


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 3-D tetrahedral cells.
class pylith::faults::TestFaultCohesiveKinTet4f : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinTet4f

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinTet4f );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinTet4f


// ----------------------------------------------------------------------
/// C++ unit testing for FaultCohesiveKin for mesh with 3-D hex cells.
class pylith::faults::TestFaultCohesiveKinHex8 : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinHex8

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestFaultCohesiveKinHex8 );

  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testIntegrateResidual );
  CPPUNIT_TEST( testIntegrateJacobian );
  CPPUNIT_TEST( testIntegrateJacobianLumped );
  CPPUNIT_TEST( testAdjustSolnLumped );
  CPPUNIT_TEST( testCalcTractionsChange );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinHex8


#endif // pylith_faults_testfaultcohesivecases_hh


// End of file 
