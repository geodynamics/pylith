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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_meshio_datawriterhdf5datamatmeshline2_hh)
#define pylith_meshio_datawriterhdf5datamatmeshline2_hh

#include "DataWriterData.hh" // ISA DataWriterData

namespace pylith {
  namespace meshio {
     class DataWriterHDF5DataMatMeshLine2;
  } // meshio
} // pylith

class pylith::meshio::DataWriterHDF5DataMatMeshLine2 : public DataWriterData
{ // DataWriterHDF5DataMatMeshLine2

public: 

  /// Constructor
  DataWriterHDF5DataMatMeshLine2(void);

  /// Destructor
  ~DataWriterHDF5DataMatMeshLine2(void);

private:

  static const char* _meshFilename; ///< Name of mesh file.
  static const char* _cellsLabel; ///< Label defining subset of cells.
  static const int _labelId; /// Value for label defining subset of cells.
  static const char* _faultLabel; ///< Name of group of vertices for fault.
  static const int _faultId; ///< Material identifier for fault.

  static const char* _timestepFilename; ///< Name of HDF5 file without fields.
  static const char* _vertexFilename; ///< Name of HDF5 file for vertex fields.
  static const char* _cellFilename; ///< Name of HDF5 file for cell fields.

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

}; // DataWriterHDF5DataMatMeshLine2

#endif // pylith_meshio_datawriterhdf5datamatmeshline2_hh

// End of file
