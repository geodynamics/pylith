// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#if !defined(pylith_meshio_datawriterdata_hh)
#define pylith_meshio_datawriterdata_hh

#include "pylith/topology/FieldBase.hh" // USES VectorFieldEnum

namespace pylith {
  namespace meshio {
     class DataWriterData;
  } // meshio
} // pylith

class pylith::meshio::DataWriterData
{ // DataWriterData

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  DataWriterData(void);

  /// Destructor
  virtual
  ~DataWriterData(void);

// PUBLIC STRUCTS ///////////////////////////////////////////////////////
public:

  struct FieldStruct {
    const char* name; ///< Name of field
    topology::FieldBase::VectorFieldEnum field_type; ///< Type of field.
    int fiber_dim; ///< Fiber dimension for field.
  }; // FieldStruct

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Name of mesh file.
  char* faultLabel; ///< Name of group of vertices for fault.
  int faultId; ///< Material identifier for fault.
  char* bcLabel; ///< Name of group of vertices for bc.

  char* timestepFilename; ///< Name of file without fields.
  char* vertexFilename; ///< Name of file for vertex fields.
  char* cellFilename; ///< Name of file for cell fields.

  PylithScalar time; ///< Time for fields.
  char* timeFormat; ///< Format for time stamp.

  char* cellsLabel; ///< Name of label for mesh cells (if using subset or boundary).
  int labelId; ///< Id for label associated with cells (if cellsLabel != 0)

  /// @name Vertex field information.
  //@{
  static const int numVertexFields; ///< Number of vertex fields.
  int numVertices; ///< Number of vertices.
  FieldStruct* vertexFieldsInfo; ///< Array of vertex field information.
  PylithScalar* vertexFields[4]; ///< Array of vertex field values.
  //@}

  /// @name Cell field information.
  //@{
  static const int numCellFields; ///< Number of cell fields.
  int numCells; ///< Number of vertices.
  FieldStruct* cellFieldsInfo; ///< Array of cell fields information.
  PylithScalar* cellFields[4]; /// Array of cell field values.
  //@}

}; // DataWriterData

#endif // pylith_meshio_datawriterdata_hh


// End of file
