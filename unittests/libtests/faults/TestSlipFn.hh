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
 * @file unittests/libtests/faults/TestSlipFn.hh
 *
 * @brief C++ TestBruneSlipFn object
 *
 * C++ unit testing for SlipFn.
 */

#if !defined(pylith_faults_testslipfn_hh)
#define pylith_faults_testslipfn_hh

#include "pylith/topology/topologyfwd.hh" // USES Mesh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestSlipFn;
  } // faults
} // pylith

/// C++ unit testing for SlipFn
class pylith::faults::TestSlipFn : public CppUnit::TestFixture
{ // class TestSlipFn

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Adjust topology of domain mesh and create fault mesh.
   *
   * @param faultMesh Fault mesh.
   * @param mesh Domain mesh.
   * @param faultLabel Label for fault.
   * @param faultId Material id for fault.
   */
  static
  void
  _createFaultMesh(topology::Mesh* faultMesh,
		   topology::Mesh* mesh,
		   const char* faultLabel,
		   const int faultId);

}; // class TestSlipFn

#endif // pylith_faults_testslipfn_hh


// End of file 
