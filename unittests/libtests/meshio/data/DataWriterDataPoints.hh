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

#if !defined(pylith_meshio_datawriterdatapoints_hh)
#define pylith_meshio_datawriterdatapoints_hh

// Data for testing writing interpolation of solution to points.

#include "DataWriterData.hh" // ISA DataWriterData

namespace pylith {
  namespace meshio {
     class DataWriterDataPoints;
  } // meshio
} // pylith

class pylith::meshio::DataWriterDataPoints : public DataWriterData
{ // DataWriterData

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  DataWriterDataPoints(void);

  /// Destructor
  virtual
  ~DataWriterDataPoints(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int numPoints; ///< Number of points for interpolation.
  int spaceDim; ///< Spatial dimension.
  PylithScalar* points; /// Points for interpolation.

}; // DataWriterData

#endif // pylith_meshio_datawriterdatapoints_hh


// End of file
