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

#if !defined(pylith_bc_DirichletPointsDatamultitet4_hh)
#define pylith_bc_DirichletPointsDatamultitet4_hh

#include "DirichletDataMulti.hh"

namespace pylith {
  namespace bc {
     class DirichletDataMultiTet4;
  } // pylith
} // bc

class pylith::bc::DirichletDataMultiTet4 : public DirichletDataMulti
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  DirichletDataMultiTet4(void);

  /// Destructor
  ~DirichletDataMultiTet4(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _numDOF; ///< Number of degrees of freedom at each point.

  static const int _numFixedDOFA; ///< Number of fixedDOF at constrained points.
  static const int _numConstrainedPtsA; ///< Number of points constrained.
  static const int _idA; ///< Boundary condition identifier
  static const char* _labelA; /// Label for boundary condition group
  static const int _fixedDOFA[]; ///< Degrees of freedom constrained at points
  static const int _constrainedPointsA[]; ///< Array of indices of constrained pts.
  static const char* _dbFilenameA; ///< Filename for db of initial values.
  static const char* _dbFilenameARate; ///< Filename for db of rate of change.
  static const PylithScalar _tRefA; ///< Reference time for rate of change.

  static const int _numFixedDOFB; ///< Number of fixedDOF at constrained points.
  static const int _numConstrainedPtsB; ///< Number of points constrained.
  static const int _idB; ///< Boundary condition identifier
  static const char* _labelB; /// Label for boundary condition group
  static const int _fixedDOFB[]; ///< Degrees of freedom constrained at points
  static const int _constrainedPointsB[]; ///< Array of indices of constrained pts.
  static const char* _dbFilenameB; ///< Filename for db of initial values.
  static const char* _dbFilenameBRate; ///< Filename for db of rate of change.
  static const PylithScalar _tRefB; ///< Reference time for rate of change.

  static const int _numFixedDOFC; ///< Number of fixedDOF at constrained points.
  static const int _numConstrainedPtsC; ///< Number of points constrained.
  static const int _idC; ///< Boundary condition identifier
  static const char* _labelC; /// Label for boundary condition group
  static const int _fixedDOFC[]; ///< Degrees of freedom constrained at points
  static const int _constrainedPointsC[]; ///< Array of indices of constrained pts.
  static const char* _dbFilenameC; ///< Filename for db of initial values.
  static const char* _dbFilenameCRate; ///< Filename for db of rate of change.
  static const PylithScalar _tRefC; ///< Reference time for rate of change.

  static const PylithScalar _field[]; ///< Values in field
  static const PylithScalar _fieldIncr[]; ///< Increment values in field
  static const int _constraintSizes[]; ///< Number of constrained DOF at each vertex
  static const int _constrainedDOF[]; ///< Indices of constrained DOF at each constrained vertex

  static const char* _meshFilename; ///< Filename of input mesh.
};

#endif // pylith_bc_dirichletdatamultitet4_hh

// End of file
