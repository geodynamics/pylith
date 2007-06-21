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

#if !defined(pylith_faults_cohesivekindata_hh)
#define pylith_faults_cohesivekindata_hh

namespace pylith {
  namespace faults {
     class CohesiveKinData;
  } // pylith
} // faults

class pylith::faults::CohesiveKinData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  CohesiveKinData(void);

  /// Destructor
  ~CohesiveKinData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Filename for input mesh

  /// @name Quadrature information
  //@{
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int cellDim; ///< Number of dimensions associated with cell
  int numBasis; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points
  double* quadPts; ///< Coordinates of quad pts in ref cell
  double* quadWts; ///< Weights of quadrature points
  double* basis; ///< Basis fns at quadrature points
  double* basisDeriv; ///< Derivatives of basis fns at quad pts
  double* verticesRef; ///< Coordinates of vertices in ref cell (dual basis)
  //@}

  /// @name Fault information
  //@{
  int id; ///< Fault material identifier
  char* label; ///< Label for fault
  char* finalSlipFilename; ///< Name of db for final slip
  char* slipTimeFilename; ///< Name of db for slip time
  char* peakRateFilename; ///< Name of db for peak rate
  char* matPropsFilename; ///< Name of db for bulk material properties
  //@}

  /// @name Input fields
  //@{
  double* fieldT; ///< Solution field at time t.
  //@}

  /// @name Calculated values.
  //@{
  double* orientation; ///< Expected values for fault orientation.
  int* constraintVertices; ///< Expected points for constraint vertices
  int* constraintCells; ///< Expected cells for constraint vertices
  double* valsResidual; ///< Expected values from residual calculation.
  double* valsJacobian; ///< Expected values from Jacobian calculation.
  double pseudoStiffness; ///< Fake stiffness for conditioning
  int numConstraintVert; ///< Number of constraint vertices
  //@}

};

#endif // pylith_faults_cohesivekindata_hh

// End of file
