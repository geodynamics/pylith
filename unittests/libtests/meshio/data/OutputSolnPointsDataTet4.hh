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

#if !defined(pylith_meshio_outputsolnpointstet4_hh)
#define pylith_meshio_outputsolnpointstet4_hh

#include "OutputSolnPointsData.hh"

namespace pylith {
  namespace meshio {
     class OutputSolnPointsDataTet4;
  } // pylith
} // meshio

class pylith::meshio::OutputSolnPointsDataTet4 : public OutputSolnPointsData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  OutputSolnPointsDataTet4(void);

  /// Destructor
  ~OutputSolnPointsDataTet4(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const char* _meshFilename; ///< Filename of input mesh

  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _numPoints; ///< Number of points.
  static const PylithScalar _points[]; ///< Coordinates of points.

  static const int _numVertices; ///< Number of points in mesh.
  static const int _fiberDim; ///< Fiber dimension of field.
  static const PylithScalar _field[]; ///< Field over mesh.
  static const PylithScalar _fieldInterp[]; ///< Field interpolated to points.

};

#endif // pylith_meshio_outputsolnpointstet4_hh


// End of file
