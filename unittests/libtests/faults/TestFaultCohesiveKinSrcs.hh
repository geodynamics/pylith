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
 * @file unittests/libtests/faults/TestFaultCohesiveKinSrcs.hh
 *
 * @brief C++ TestFaultCohesiveKinSrcs object
 *
 * C++ unit testing for FaultCohesiveKin.
 */

#if !defined(pylith_faults_testfaultcohesivekinsrcs_hh)
#define pylith_faults_testfaultcohesivekinsrcs_hh

#include "TestFaultCohesiveKin.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveKinSrcs;
  } // faults
} // pylith

/// C++ unit testing for FaultCohesiveKin
class pylith::faults::TestFaultCohesiveKinSrcs : public TestFaultCohesiveKin
{ // class TestFaultCohesiveKinSrcs

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Setup testing data.
  void setUp(void);

}; // class TestFaultCohesiveKinSrcs

#endif // pylith_faults_testfaultcohesivekinsrcs_hh


// End of file 
