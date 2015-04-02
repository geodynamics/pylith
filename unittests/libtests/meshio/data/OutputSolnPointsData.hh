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

#if !defined(pylith_meshio_outputsolnpointsdata_hh)
#define pylith_meshio_outputsolnpointsdata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace meshio {
     class OutputSolnPointsData;
  } // pylith
} // meshio

class pylith::meshio::OutputSolnPointsData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  OutputSolnPointsData(void);

  /// Destructor
  ~OutputSolnPointsData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Filename for input mesh

  /// @name Point information.
  //@{
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int numPoints; ///< Number of points.
  PylithScalar* points; ///< Coordinates of points.
  //@}

  /// @name Input data.
  //@{
  int numVertices;
  int fiberDim;
  PylithScalar* field; ///< Field over mesh.
  //@}

  /// @name Calculated values.
  //@{
  PylithScalar* fieldInterp; ///< Field interpolate to points.
  //@}

};

#endif // pylith_meshio_outputsolnpointsdata_hh

// End of file
