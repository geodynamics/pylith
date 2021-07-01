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

#if !defined(pylith_meshio_outputsolnpointshex8_hh)
#define pylith_meshio_outputsolnpointshex8_hh

#include "OutputSolnPointsData.hh"

namespace pylith {
  namespace meshio {
     class OutputSolnPointsDataHex8;
  } // pylith
} // meshio

class pylith::meshio::OutputSolnPointsDataHex8 : public OutputSolnPointsData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  OutputSolnPointsDataHex8(void);

  /// Destructor
  ~OutputSolnPointsDataHex8(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const char* _meshFilename; ///< Filename of input mesh

    static const int _spaceDim; ///< Number of dimensions in vertex coordinates
    static const int _numPoints; ///< Number of points.
    static const PylithScalar _points[]; ///< Coordinates of points.
    static const char* _names[]; /// Names of stations.

    static const int _fiberDim; ///< Fiber dimension of field.
    static const PylithScalar _coefs[]; ///< Polynomial coefficients for field.

};

#endif // pylith_meshio_outputsolnpointshex8_hh


// End of file
