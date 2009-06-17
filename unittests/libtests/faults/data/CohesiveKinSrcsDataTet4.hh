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

#if !defined(pylith_faults_cohesivekindatatet4_hh)
#define pylith_faults_cohesivekindatatet4_hh

#include "CohesiveKinData.hh"

namespace pylith {
  namespace faults {
     class CohesiveKinSrcsDataTet4;
  } // pylith
} // faults

class pylith::faults::CohesiveKinSrcsDataTet4 : public CohesiveKinData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  CohesiveKinSrcsDataTet4(void);

  /// Destructor
  ~CohesiveKinSrcsDataTet4(void);

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

  static const double _fieldT[]; ///< Solution field at time t.

  static const double _orientation[]; ///< Expected values for fault orientation.
  static const double _area[]; ///< Expected values for fault area.
  static const int _constraintVertices[]; ///< Expected points for constraint vertices
  static const int _constraintCells[]; ///< Expected cells for constraint vertices
  static const double _valsResidual[]; ///< Expected values from residual calculation.
  static const double _valsResidualIncr[]; ///< Expected values from residual calculation using solution increment.
  static const double _valsJacobian[]; ///< Expected values from Jacobian calculation.
  static const int _numConstraintVert; ///< Number of constraint vertices

};

#endif // pylith_faults_cohesivekindatatet4_hh


// End of file
