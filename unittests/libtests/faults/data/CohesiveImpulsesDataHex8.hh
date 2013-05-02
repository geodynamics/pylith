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

#if !defined(pylith_faults_cohesiveimpulsesdatahex8_hh)
#define pylith_faults_cohesiveimpulsesdatahex8_hh

#include "CohesiveImpulsesData.hh"

namespace pylith {
  namespace faults {
     class CohesiveImpulsesDataHex8;
  } // pylith
} // faults

class pylith::faults::CohesiveImpulsesDataHex8 : public CohesiveImpulsesData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  CohesiveImpulsesDataHex8(void);

  /// Destructor
  ~CohesiveImpulsesDataHex8(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const char* _meshFilename; ///< Filename of input mesh

  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _cellDim; ///< Number of dimensions associated with cell

  static const int _numBasis; ///< Number of vertices in cell
  static const int _numQuadPts; ///< Number of quadrature points
  static const PylithScalar _quadPts[]; ///< Coordinates of quad pts in ref cell
  static const PylithScalar _quadWts[]; ///< Weights of quadrature points
  static const PylithScalar _basis[]; ///< Basis fns at quadrature points
  static const PylithScalar _basisDeriv[]; ///< Derivatives of basis fns at quad pts
  static const PylithScalar _verticesRef[]; ///< Coordinates of vertices in ref cell (dual basis)

  static const int _id; ///< Fault material identifier
  static const char* _label; ///< Label for fault
  static const char* _impulseAmpFilename; ///< Name of db for impulse amplitude
  static const int _impulseDOF[]; ///< Fault components associated with impulses.
  static const int _numComponents; ///< Number of components;
  //@}

  static const PylithScalar _fieldT[]; ///< Field over domain at time t.
  static const PylithScalar _fieldIncr[]; ///< Solution increment field over domain at time t.

  static const PylithScalar _orientation[]; ///< Expected values for fault orientation.
  static const PylithScalar _area[]; ///< Expected values for fault area.
  static const PylithScalar _amplitude[]; ///< Expected values for impulse amplitude.
  static const int _numImpulses; ///< Number of impulses.
  static const PylithScalar _residual[]; ///< Expected values from residual.

  static const int _constraintVertices[]; ///< Expected points for constraint vertices
  static const int _negativeVertices[]; ///< Expected points for negative-side fault vertices
  static const int _numConstraintVert; ///< Number of constraint vertices

};

#endif // pylith_faults_cohesiveimpulsesdatahex8_hh


// End of file
