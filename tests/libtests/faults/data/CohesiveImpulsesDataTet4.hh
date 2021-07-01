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

#if !defined(pylith_faults_cohesiveimpulsesdatatet4_hh)
#define pylith_faults_cohesiveimpulsesdatatet4_hh

#include "CohesiveImpulsesData.hh"

namespace pylith {
  namespace faults {
     class CohesiveImpulsesDataTet4;
  } // pylith
} // faults

class pylith::faults::CohesiveImpulsesDataTet4 : public CohesiveImpulsesData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  CohesiveImpulsesDataTet4(void);

  /// Destructor
  ~CohesiveImpulsesDataTet4(void);

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
  static const PylithScalar _residual[]; ///< Expected values from residual calculation.

  static const int _constraintEdges[]; ///< Expected points for constraint edges
  static const int _negativeVertices[]; ///< Expected points for negative-side fault vertices
  static const int _numConstraintEdges; ///< Number of constraint edges

};

#endif // pylith_faults_cohesiveimpulsesdatatet4_hh


// End of file
