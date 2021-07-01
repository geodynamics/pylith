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

#if !defined(pylith_bc_dirichletdatahex8_hh)
#define pylith_bc_dirichletdatahex8_hh

#include "DirichletData.hh"

namespace pylith {
  namespace bc {
     class DirichletDataHex8;
  } // pylith
} // bc

class pylith::bc::DirichletDataHex8 : public DirichletData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  DirichletDataHex8(void);

  /// Destructor
  ~DirichletDataHex8(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _numDOF; ///< Number of degrees of freedom at each point.

  static const int _numFixedDOF; ///< Number of fixedDOF at constrained points.
  static const int _numConstrainedPts; ///< Number of points constrained.

  static const int _id; ///< Boundary condition identifier
  static const char* _label; /// Label for boundary condition group

  static const int _fixedDOF[]; ///< Degrees of freedom constrained at points

  static const int _constrainedPoints[]; ///< Array of indices of constrained pts.
  static const PylithScalar _tRef; ///< Reference time for rate of change of value
  static const PylithScalar _valuesInitial[]; ///< Initial values.
  static const PylithScalar _valueRate; ///< Rate of change of values.

  static const char* _meshFilename; ///< Filename of input mesh.
  static const char* _dbFilename; ///< Filename of simple spatial database.
};

#endif // pylith_bc_dirichletdatahex8_hh

// End of file
