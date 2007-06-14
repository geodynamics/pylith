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

#if !defined(pylith_feassemble_integratordata_hh)
#define pylith_feassemble_integratordata_hh

namespace pylith {
  namespace feassemble {
     class IntegratorData;
  } // pylith
} // feassemble

class pylith::feassemble::IntegratorData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  IntegratorData(void);

  /// Destructor
  ~IntegratorData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  /// @name Mesh information
  //@{
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int cellDim; ///< Number of dimensions associated with cell
  int numVertices; ///< Number of vertices
  int numCells; ///< Number of cells
  double* vertices; ///< Coordinates of vertices
  int* cells; ///< Indices of vertices in cells
  double* verticesRef; ///< Coordinates of vertices in ref cell (dual basis)
  //@}

  /// @name Quadrature information
  //@{
  int numBasis; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points
  double* quadPts; ///< Coordinates of quad pts in ref cell
  double* quadWts; ///< Weights of quadrature points
  double* basis; ///< Basis fns at quadrature points
  double* basisDerivRef; ///< Derivatives of basis fns at quad pts
  //@}

  /// @name Material information
  //@{
  char* matType; ///< String corresponding to material type.
  char* matDBFilename; ///< Filename for database of material properties.
  int matId; ///< Material identifier.
  char* matLabel; ///< Label of material.
  //@}

  /// @name Input fields
  //@{
  double dt; ///< Time step
  double* fieldTpdt; ///< Input field at time t+dt.
  double* fieldT; ///< Input field at time t.
  double* fieldTmdt; ///< Input field at time t-dt.
  //@}

  /// @name Calculated values.
  //@{
  double* valsResidual; ///< Expected values from residual calculation.
  double* valsJacobian; ///< Expected values from Jacobian calculation.
  //@}
};

#endif // pylith_feassemble_integratordata_hh

// End of file
