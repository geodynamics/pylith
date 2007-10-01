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

#if !defined(pylith_bc_absorbingdampersdata_hh)
#define pylith_bc_absorbingdampersdata_hh

namespace pylith {
  namespace bc {
     class AbsorbingDampersData;
  } // pylith
} // bc

class pylith::bc::AbsorbingDampersData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  AbsorbingDampersData(void);

  /// Destructor
  ~AbsorbingDampersData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Name of file with input mesh

  /// @name Boundary mesh information
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

  /// @name Parameter information
  //@{
  char* spatialDBFilename; ///< Filename for database of parameters.
  int id; ///< Identifier for boundary condition
  char* label; ///< Label for boundary condition
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

#endif // pylith_bc_absorbingdampersdata_hh


// End of file
