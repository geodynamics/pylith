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
 * @file unittests/libtests/faults/TestFaultMesh.hh
 *
 * @brief C++ object for constructing fault mesh.
 */

#if !defined(pylith_faults_testfaultmesh_hh)
#define pylith_faults_testfaultmesh_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultMesh;
  } // faults
} // pylith

#include "pylith/topology/topologyfwd.hh" // USES Mesh

/// C++ object for constructing fault mesh.
class pylith::faults::TestFaultMesh
{ // class TestFaultMesh

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Adjust topology of domain mesh and create fault mesh.
   *
   * @param faultMesh Fault mesh.
   * @param mesh Domain mesh.
   * @param faultLabel Label for fault.
   * @param faultId Material id for fault.
   */
  static
  void createFaultMesh(topology::Mesh* faultMesh,
		       topology::Mesh* mesh,
		       const char* faultLabel,
		       const int faultId);

}; // class TestFaultMesh

#endif // pylith_faults_testfaultmesh_hh

// End of file 
