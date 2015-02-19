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
 * @file unittests/libtests/meshio/TestDataWriterPoints.hh
 *
 * @brief C++ TestDataWriterPoints object
 *
 * C++ unit testing for DataWriter<Mesh>.
 */

#if !defined(pylith_meshio_testdatawriterpoints_hh)
#define pylith_meshio_testdatawriterpoints_hh

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterPoints;

    class DataWriterDataPoints;
  } // meshio
} // pylith

/// C++ unit testing for iDataWriter<Mesh> with interpolated solution.
class pylith::meshio::TestDataWriterPoints
{ // class TestDataWriterPoints

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Initialize mesh.
  void _initialize(void);

  /** Create vertex fields.
   *
   * @param fields Vertex fields.
   */
  void
  _createVertexFields(topology::Fields* fields) const;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DataWriterDataPoints* _data; ///< Data for testing
  topology::Mesh* _mesh; ///< Mesh for data

}; // class TestDataWriterPoints

#endif // pylith_meshio_testdatawriterpoints_hh


// End of file 
