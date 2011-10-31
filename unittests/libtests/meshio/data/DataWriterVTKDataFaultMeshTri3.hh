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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_meshio_datawritervtkdatafaultmeshtri3_hh)
#define pylith_meshio_datawritervtkdatafaultmeshtri3_hh

#include "DataWriterData.hh" // ISA DataWriterData

namespace pylith {
  namespace meshio {
     class DataWriterVTKDataFaultMeshTri3;
  } // meshio
} // pylith

class pylith::meshio::DataWriterVTKDataFaultMeshTri3 : public DataWriterData
{ // DataWriterVTKDataFaultMeshTri3

public: 

  /// Constructor
  DataWriterVTKDataFaultMeshTri3(void);

  /// Destructor
  ~DataWriterVTKDataFaultMeshTri3(void);

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

  static const PylithScalar _vertexField0[]; ///< Values for vertex field 0.
  static const PylithScalar _vertexField1[]; ///< Values for vertex field 1.
  static const PylithScalar _vertexField2[]; ///< Values for vertex field 2.
  //@}

  /// @name Cell field information.
  //@{
  static const int _numCellFields; ///< Number of cell fields.
  static const int _numCells; ///< Number of cells.
  static const FieldStruct _cellFields[]; ///< Array of cell fields.

  static const PylithScalar _cellField0[]; ///< Values for cell field 0.
  static const PylithScalar _cellField1[]; ///< Values for cell field 1.
  static const PylithScalar _cellField2[]; ///< Values for cell field 2.
  //@}

}; // DataWriterVTKDataFaultMeshTri3

#endif // pylith_meshio_datawritervtkdatafaultmeshtri3_hh

// End of file
