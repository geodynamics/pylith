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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_bc_forcedatadatatri3_hh)
#define pylith_bc_forcedatadatatri3_hh

#include "PointForceData.hh"

namespace pylith {
  namespace bc {
     class PointForceDataTri3;
  } // pylith
} // bc

class pylith::bc::PointForceDataTri3 : public PointForceData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  PointForceDataTri3(void);

  /// Destructor
  ~PointForceDataTri3(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const double _tRef; ///< Reference time for rate of change of forces.
  static const double _forceRate; ///< Rate of change of force.
  static const double _tResidual; ///< Time for computing residual.

  static const int _numDOF; ///< Number of degrees of freedom at each point.
  static const int _numForceDOF; ///< Number of forces at points.
  static const int _numForcePts; ///< Number of points with forces.

  static const int _id; ///< Boundary condition identifier
  static const char* _label; ///< Label for boundary condition group

  static const int _forceDOF[]; ///< Degrees of freedom that are constrained at each point
  static const int _forcePoints[]; ///< Array of indices of points with forces.
  static const double _forceInitial[]; ///< Forces at points.
  static const double _residual[]; ///< Residual field.

  static const char* _meshFilename; ///< Filename for input mesh.
  static const char* _dbFilename; ///< Filename of simple spatial database.
};

#endif // pylith_bc_forcedatadatatri3_hh

// End of file
