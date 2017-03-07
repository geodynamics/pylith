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

#if !defined(pylith_materials_testsubmesh_data_hh)
#define pylith_materials_testsubmesh_data_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
    namespace topology {
	class TestSubMesh_Data;
    } // topology
} // pylith

class pylith::topology::TestSubMesh_Data
{ // TestSubMesh_Data
    
// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  TestSubMesh_Data(void);

  /// Destructor
  ~TestSubMesh_Data(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // GENERAL, VALUES DEPEND ON TEST CASE

    /// @defgroup Domain mesh information.
    /// @{
    int cellDim; ///< Cell dimension (matches space dimension).
    int numVertices; ///< Number of vertices.
    int numCells; ///< Number of cells.
    int numCorners; ///< Number of vertices per cell.
    int* cells; ///< Array of vertices in cells [numCells*numCorners].
    PylithScalar* coordinates; ///< Coordinates of vertices [numVertices*cellDim].
    const char* label; ///< Label of group associated with submesh.
    int groupSize; ///< Number of vertices in group.
    int* groupVertices; ///< Array of vertices in group.
    /// @}

    /// @defgroup SubMesh information.
    /// @{
    int submeshNumCorners; ///< Number of vertices per cell.
    int submeshNumVertices; ///< Number of vertices in submesh.
    int* submeshVertices; ///< Vertices in submesh.
    int submeshNumCells; ///< Number of cells in submesh.
    int* submeshCells; ///< Array of vertices in cells [submeshNumCells*submeshNumCorners].
    /// @}
  
}; // TestSubMesh_Data

#endif // pylith_topology_testsubmesh_data_hh

// End of file
