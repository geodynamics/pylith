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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#if !defined(pylith_faults_cohesivekindatahex8_hh)
#define pylith_faults_cohesivekindatahex8_hh

#include "CohesiveKinData.hh"

namespace pylith {
  namespace faults {
     class CohesiveKinDataHex8;
  } // pylith
} // faults

class pylith::faults::CohesiveKinDataHex8 : public CohesiveKinData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  CohesiveKinDataHex8(void);

  /// Destructor
  ~CohesiveKinDataHex8(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const char* _meshFilename; ///< Filename of input mesh

  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _cellDim; ///< Number of dimensions associated with cell

  static const int _numBasis; ///< Number of vertices in cell
  static const int _numQuadPts; ///< Number of quadrature points
  static const double _quadPts[]; ///< Coordinates of quad pts in ref cell
  static const double _quadWts[]; ///< Weights of quadrature points
  static const double _basis[]; ///< Basis fns at quadrature points
  static const double _basisDeriv[]; ///< Derivatives of basis fns at quad pts
  static const double _verticesRef[]; ///< Coordinates of vertices in ref cell (dual basis)

  static const int _id; ///< Fault material identifier
  static const char* _label; ///< Label for fault
  static const char* _finalSlipFilename; ///< Name of db for final slip
  static const char* _slipTimeFilename; ///< Name of db for slip time
  static const char* _riseTimeFilename; ///< Name of db for rise time
  static const char* _matPropsFilename; ///< Name of db for bulk mat properties.
  //@}

  static const double _fieldT[]; ///< Field over domain at time t.
  static const double _fieldIncr[]; ///< Solution increment field over domain at time t.
  static const double _jacobianLumped[]; ///< Lumped Jacobian.

  static const double _orientation[]; ///< Expected values for fault orientation.
  static const double _area[]; ///< Expected values for fault area.
  static const double _residual[]; ///< Expected values from residual calculation.
  static const double _residualIncr[]; ///< Expected values from residual calculation with solution increment.
  static const double _jacobian[]; ///< Expected values from Jacobian calculation.
  static const double _fieldIncrAdjusted[]; ///< Expected values for colution increment after adjustment.

  static const int _verticesFault[]; ///< Expected points for Fault vertices
  static const int _verticesLagrange[]; ///< Expected points for Lagrange vertices
  static const int _verticesPositive[]; ///< Expected points for vertices on + side of fault.
  static const int _verticesNegative[]; ///< Expected points for vertices on - side of fault.
  static const int _numFaultVertices; ///< Number of fault vertices

  static const int _numCohesiveCells; ///< Number of cohesive cells
  static const int _cellMappingFault[]; ///< Fault cell
  static const int _cellMappingCohesive[]; ///< Cohesive cell

};

#endif // pylith_faults_cohesivekindatahex8_hh


// End of file
