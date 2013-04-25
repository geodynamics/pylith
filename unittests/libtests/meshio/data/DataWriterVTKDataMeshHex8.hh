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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_meshio_datawritervtkdatameshhex8_hh)
#define pylith_meshio_datawritervtkdatameshhex8_hh

#include "DataWriterData.hh" // ISA DataWriterData

namespace pylith {
  namespace meshio {
     class DataWriterVTKDataMeshHex8;
  } // meshio
} // pylith

class pylith::meshio::DataWriterVTKDataMeshHex8 : public DataWriterData
{ // DataWriterVTKDataMeshHex8

public: 

  /// Constructor
  DataWriterVTKDataMeshHex8(void);

  /// Destructor
  ~DataWriterVTKDataMeshHex8(void);

private:

  static const char* _meshFilename; ///< Name of mesh file.
  static const char* _faultLabel; ///< Name of group of vertices for fault.
  static const int _faultId; ///< Material identifier for fault.

  static const char* _timestepFilename; ///< Name of VTK file without fields.
  static const char* _vertexFilename; ///< Name of VTK file for vertex fields.
  static const char* _cellFilename; ///< Name of VTK file for cell fields.

  static const PylithScalar _time; ///< Time for fields.
  static const char* _timeFormat; ///< Format for time stamp.

  /// @name Vertex field information.
  //@{
  static const int _numVertexFields; ///< Number of vertex fields.
  static const int _numVertices; ///< Number of vertices.
  static const FieldStruct _vertexFields[]; ///< Array of vertex fields.

  static const PylithScalar _vertexFieldScalar[]; ///< Values for scalar vertex field.
  static const PylithScalar _vertexFieldVector[]; ///< Values for vector vertex field .
  static const PylithScalar _vertexFieldTensor[]; ///< Values for tensor vertex field.
  static const PylithScalar _vertexFieldOther[]; ///< Values for other vertex field.
  //@}

  /// @name Cell field information.
  //@{
  static const int _numCellFields; ///< Number of cell fields.
  static const int _numCells; ///< Number of cells.
  static const FieldStruct _cellFields[]; ///< Array of cell fields.

  static const PylithScalar _cellFieldScalar[]; ///< Values for scalar cell field.
  static const PylithScalar _cellFieldVector[]; ///< Values for vector cell field.
  static const PylithScalar _cellFieldTensor[]; ///< Values for tensor cell field.
  static const PylithScalar _cellFieldOther[]; ///< Values for other cell field.
  //@}

}; // DataWriterVTKDataMeshHex8

#endif // pylith_meshio_datawritervtkdatameshhex8_hh

// End of file
