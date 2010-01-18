// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#if !defined(pylith_faults_cohesivedynldataquad4_hh)
#define pylith_faults_cohesivedynldataquad4_hh

#include "CohesiveDynLData.hh"

namespace pylith {
  namespace faults {
     class CohesiveDynLDataQuad4;
  } // pylith
} // faults

class pylith::faults::CohesiveDynLDataQuad4 : public CohesiveDynLData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  CohesiveDynLDataQuad4(void);

  /// Destructor
  ~CohesiveDynLDataQuad4(void);

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
  static const char* _initialTractFilename; ///< Name of db for initial tractions.

  static const double _fieldT[]; ///< Solution field at time t.
  static const double _fieldIncrStick[]; ///< Solution increment at time t for stick case.
  static const double _fieldIncrSlip[]; ///< Solution increment at time t for slip case.
  static const double _fieldIncrOpen[]; ///< Solution increment at time t for opening case.
  static const double _jacobian[]; ///< Jacobian sparse matrix.

  static const double _orientation[]; ///< Expected values for fault orientation.
  static const double _area[]; ///< Expected values for fault area.
  static const double _initialTractions[]; ///< Expected values for initial tractions.
  static const double _fieldIncrSlipE[]; ///< Expected values for solution increment for slip case.
  static const double _slipSlipE[]; ///< Expected values for slip for slip case.
  static const double _fieldIncrOpenE[]; ///< Expected values for solution increment for opening case.
  static const double _slipOpenE[]; ///< Expected values for slip for opening case.
  static const int _constraintVertices[]; ///< Expected points for constraint vertices
  static const int _constraintCells[]; ///< Expected cells for constraint vertices
  static const int _numConstraintVert; ///< Number of constraint vertices

};

#endif // pylith_faults_cohesivedynldataquad4_hh


// End of file
