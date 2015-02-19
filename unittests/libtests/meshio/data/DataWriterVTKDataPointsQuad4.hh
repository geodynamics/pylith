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

#if !defined(pylith_meshio_datawritervtkdatapointsquad4_hh)
#define pylith_meshio_datawritervtkdatapointsquad4_hh

#include "DataWriterDataPoints.hh" // ISA DataWriterData

namespace pylith {
  namespace meshio {
     class DataWriterVTKDataPointsQuad4;
  } // meshio
} // pylith

class pylith::meshio::DataWriterVTKDataPointsQuad4 : public DataWriterDataPoints
{ // DataWriterVTKDataPointsQuad4

public: 

  /// Constructor
  DataWriterVTKDataPointsQuad4(void);

  /// Destructor
  ~DataWriterVTKDataPointsQuad4(void);

private:

  static const char* _meshFilename; ///< Name of mesh file.

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

}; // DataWriterVTKDataPointsQuad4

#endif // pylith_meshio_datawritervtkdatapointsquad4_hh

// End of file
