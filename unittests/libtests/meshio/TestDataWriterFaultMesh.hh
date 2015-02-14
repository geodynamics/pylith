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
 * @file unittests/libtests/meshio/TestDataWriterFaultMesh.hh
 *
 * @brief C++ TestDataWriterFaultMesh object
 *
 * C++ unit testing for DataWriter<FaultMesh>.
 */

#if !defined(pylith_meshio_testdatawriterfaultmesh_hh)
#define pylith_meshio_testdatawriterfaultmesh_hh

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterFaultMesh;

    class DataWriterData;
  } // meshio
} // pylith

/// C++ unit testing for DataWriter<FaultMesh>.
class pylith::meshio::TestDataWriterFaultMesh
{ // class TestDataWriterFaultMesh

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

  /** Create cell fields.
   *
   * @param fields Cell fields.
   */
  void
  _createCellFields(topology::Fields* fields) const;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DataWriterData* _data; ///< Data for testing
  topology::Mesh* _mesh; ///< Mesh for domain
  topology::Mesh* _faultMesh; ///< Fault mesh.

}; // class TestDataWriterFaultMesh

#endif // pylith_meshio_testdatawriterfaultmesh_hh


// End of file 
