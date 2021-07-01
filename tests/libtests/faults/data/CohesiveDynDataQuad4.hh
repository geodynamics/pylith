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

#if !defined(pylith_faults_cohesivedynldataquad4_hh)
#define pylith_faults_cohesivedynldataquad4_hh

#include "CohesiveDynData.hh"

namespace pylith {
  namespace faults {
     class CohesiveDynDataQuad4;
  } // pylith
} // faults

class pylith::faults::CohesiveDynDataQuad4 : public CohesiveDynData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  CohesiveDynDataQuad4(void);

  /// Destructor
  ~CohesiveDynDataQuad4(void);

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
  static const char* _initialTractFilename; ///< Name of db for initial tractions.

  static const PylithScalar _fieldT[]; ///< Solution field at time t.
  static const PylithScalar _fieldIncrStick[]; ///< Solution increment at time t for stick case.
  static const PylithScalar _fieldIncrSlip[]; ///< Solution increment at time t for slip case.
  static const PylithScalar _fieldIncrOpen[]; ///< Solution increment at time t for opening case.
  static const PylithScalar _jacobian[]; ///< Jacobian sparse matrix.

  static const PylithScalar _orientation[]; ///< Expected values for fault orientation.
  static const PylithScalar _area[]; ///< Expected values for fault area.
  static const PylithScalar _initialTractions[]; ///< Expected values for initial tractions.
  static const PylithScalar _slipStickE[]; ///< Expected values for slip for stick case.
  static const PylithScalar _fieldIncrSlipE[]; ///< Expected values for solution increment for slip case.
  static const PylithScalar _slipSlipE[]; ///< Expected values for slip for slip case.
  static const PylithScalar _fieldIncrOpenE[]; ///< Expected values for solution increment for opening case.
  static const PylithScalar _slipOpenE[]; ///< Expected values for slip for opening case.
  static const int _constraintEdges[]; ///< Expected points for constraint edges
  static const int _negativeVertices[]; ///< Expected points for negative-side fault vertices
  static const int _numConstraintEdges; ///< Number of constraint edges

};

#endif // pylith_faults_cohesivedynldataquad4_hh


// End of file
