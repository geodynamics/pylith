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

#if !defined(pylith_meshio_datawritervtkdatapointshex8_hh)
#define pylith_meshio_datawritervtkdatapointshex8_hh

#include "DataWriterDataPoints.hh" // ISA DataWriterData

namespace pylith {
  namespace meshio {
     class DataWriterVTKDataPointsHex8;
  } // meshio
} // pylith

class pylith::meshio::DataWriterVTKDataPointsHex8 : public DataWriterDataPoints
{ // DataWriterVTKDataPointsHex8

public: 

  /// Constructor
  DataWriterVTKDataPointsHex8(void);

  /// Destructor
  ~DataWriterVTKDataPointsHex8(void);

private:

  static const char* _meshFilename; ///< Name of mesh file.
  static const char* _faultLabel; ///< Name of group of vertices for fault.
  static const int _faultId; ///< Material identifier for fault.

  static const char* _timestepFilename; ///< Name of VTK file without fields.
  static const char* _vertexFilename; ///< Name of VTK file for vertex fields.

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

  /// @name Point information.
  //@{
  static const int _numPoints; ///< Number of points.
  static const int _spaceDim; ///< Spatial dimension.
  static const PylithScalar _points[]; ///< Coordinates of points.
  //@}

}; // DataWriterVTKDataPointsHex8

#endif // pylith_meshio_datawritervtkdatapointshex8_hh

// End of file
