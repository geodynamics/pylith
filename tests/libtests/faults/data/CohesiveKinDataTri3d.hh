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

#if !defined(pylith_faults_cohesivekindatatri3d_hh)
#define pylith_faults_cohesivekindatatri3d_hh

#include "CohesiveKinData.hh"

namespace pylith {
  namespace faults {
     class CohesiveKinDataTri3d;
  } // pylith
} // faults

class pylith::faults::CohesiveKinDataTri3d : public CohesiveKinData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  CohesiveKinDataTri3d(void);

  /// Destructor
  ~CohesiveKinDataTri3d(void);

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
  static const char* _finalSlipFilename; ///< Name of db for final slip
  static const char* _slipTimeFilename; ///< Name of db for slip time
  static const char* _riseTimeFilename; ///< Name of db for rise time
  static const char* _matPropsFilename; ///< Name of db for bulk mat properties.
  //@}

  static const PylithScalar _fieldT[]; ///< Solution field at time t.
  static const PylithScalar _fieldIncr[]; ///< Solution increment field over domain at time t.
  static const PylithScalar _jacobianLumped[]; ///< Lumped Jacobian.

  static const PylithScalar _orientation[]; ///< Expected values for fault orientation.
  static const PylithScalar _area[]; ///< Expected values for fault area.
  static const PylithScalar _residual[]; ///< Expected values from residual calculation.
  static const PylithScalar _jacobian[]; ///< Expected values from Jacobian calculation.

  static const int _verticesFault[]; ///< Expected points for Fault vertices
  static const int _edgesLagrange[]; ///< Expected points for Lagrange multipliers
  static const int _verticesPositive[]; ///< Expected points for vertices on + side of fault.
  static const int _verticesNegative[]; ///< Expected points for vertices on - side of fault.
  static const int _numFaultVertices; ///< Number of fault vertices

  static const int _numCohesiveCells; ///< Number of cohesive cells
  static const int _cellMappingFault[]; ///< Fault cell
  static const int _cellMappingCohesive[]; ///< Cohesive cell

};

#endif // pylith_faults_cohesivekindatatri3d_hh


// End of file
