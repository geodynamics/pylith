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

/**
 * @file unittests/libtests/meshio/TestDataWriterSubMesh.hh
 *
 * @brief C++ TestDataWriterSubMesh object
 *
 * C++ unit testing for DataWriter<SubMesh>.
 */

#if !defined(pylith_meshio_testdatawritersubmesh_hh)
#define pylith_meshio_testdatawritersubmesh_hh

#include "pylith/topology/topologyfwd.hh" // USES Mesh, SubMesh, Field

/// Namespace for pylith package
namespace pylith {
  namespace meshio {
    class TestDataWriterSubMesh;

    class DataWriterData;
  } // meshio
} // pylith

/// C++ unit testing for DataWriter<SubMesh>.
class pylith::meshio::TestDataWriterSubMesh
{ // class TestDataWriterSubMesh

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
  _createVertexFields(topology::Fields<topology::Field<topology::Mesh> >* fields) const;

  /** Create cell fields.
   *
   * @param fields Cell fields.
   */
  void
  _createCellFields(topology::Fields<topology::Field<topology::SubMesh> >* fields) const;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  DataWriterData* _data; ///< Data for testing
  topology::Mesh* _mesh; ///< Mesh for domain
  topology::SubMesh* _submesh; ///< Mesh for subdomain.
  bool _flipFault; ///< If true, flip fault orientation.

}; // class TestDataWriterSubMesh

#endif // pylith_meshio_testdatawritersubmesh_hh


// End of file 
