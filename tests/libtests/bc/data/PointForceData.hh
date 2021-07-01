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

#if !defined(pylith_bc_pointforcedata_hh)
#define pylith_bc_pointforcedata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
  namespace bc {
     class PointForceData;
  } // pylith
} // bc

class pylith::bc::PointForceData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  PointForceData(void);

  /// Destructor
  ~PointForceData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  PylithScalar tRef; ///< Reference time for rate of change of forces.
  PylithScalar forceRate; ///< Rate of change of force.
  PylithScalar tResidual; ///< Time for computing residual.

  int numDOF; ///< Number of degrees of freedom at each point.
  int numForceDOF; ///< Number of forces at points.
  int numForcePts; ///< Number of points with forces.

  int id; ///< Boundary condition identifier
  char* label; ///< Label for boundary condition group

  int* forceDOF; ///< Degrees of freedom that are constrained at each point
  int* forcePoints; ///< Array of indices of points with forces.
  PylithScalar* forceInitial; ///< Forces at points.
  PylithScalar* residual; ///< Residual field.

  char* meshFilename; ///< Filename for input mesh.
  char* dbFilename; ///< Filename of simple spatial database.

  /// @name Scales information for nondimensionalization.
  //@{
  PylithScalar lengthScale; ///< Length scale.
  PylithScalar pressureScale; ///< Pressure scale.
  PylithScalar timeScale; ///< Time scale.
  PylithScalar densityScale; ///< Density scale.
  //@}

};

#endif // pylith_bc_pointforcedata_hh

// End of file
