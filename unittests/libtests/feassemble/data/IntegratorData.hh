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

  int numVertices; ///< Number of vertices
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int numCells; ///< Number of cells
  int cellDim; ///< Number of dimensions associated with cell
  int numCorners; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points
  int fiberDim; ///< Number of values per vertex in field

  /// @name Mesh information
  //@{
  double* vertices; ///< Pointer to coordinates of vertices
  int* cells; ///< Pointer to indices of vertices in cells
  //@}

  /// @name Discretization information
  //@{
  double* quadPts; ///< Coordinates of quad pts in ref cell
  double* quadWts; ///< Weights of quadrature points
  double* basis; ///< Basis fns at quadrature points
  double* basisDeriv; ///< Derivatives of basis fns at quad pts
  //@}

  /// @name Integration information
  //@{
  double* fieldIn; ///< Input field for integration action
  double* valsAction; ///< Expected output for integration action
  double* valsMatrix; ///< Expected output for integration
  //@}
};

#endif // pylith_feassemble_integratordata_hh

// End of file
