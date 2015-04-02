// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file playpen/closure/TestClosure.hh
 *
 */

#if !defined(pylith_playpen_testclosure_hh)
#define pylith_playpen_testclosure_hh

// Include directives ---------------------------------------------------
namespace pylith {
  namespace playpen {
    class TestClosure;
  } // playpen
} // pylith

#include "pylith/topology/topologyfwd.hh"

// ElasticityExplicitTet4 ---------------------------------------------------
class pylith::playpen::TestClosure
{ // TestClosure

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  TestClosure(void);

  /// Destructor
  ~TestClosure(void);

  /** Set number of iterations.
   *
   * @param value Number of iterations.
   */
  void iterations(const long value);

  /** Test restrictClosure().
   *
   * @param mesh Finite-element mesh.
   */
  void testRestrictClosure(const pylith::topology::Mesh& mesh);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  long _niterations; ///< Number of iterations.

  static const int _spaceDim; ///< Space dimension.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  TestClosure(const TestClosure&); ///< Not implemented
  const TestClosure& operator=(const TestClosure&); ///< Not implemented.

}; // TestClosure

#endif // pylith_playpen_testclosure_hh


// End of file 
