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

#if !defined(pylith_feassemble_quadraturedata_hh)
#define pylith_feassemble_quadraturedata_hh

namespace pylith {
  namespace feassemble {
     class QuadratureData;
  } // pylith
} // feassemble

class pylith::feassemble::QuadratureData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  QuadratureData(void);

  /// Destructor
  ~QuadratureData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int numVertices; ///< Number of vertices
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int numCells; ///< Number of cells (=1)
  int cellDim; ///< Number of dimensions associated with cell
  int numBasis; ///< Number of vertices in cell
  int numQuadPts; ///< Number of quadrature points

  double* vertices; ///< Pointer to coordinates of vertices
  int* cells; ///< Pointer to indices of vertices in cells

  double* verticesRef; ///< Coordinates of vertices in ref cell
  double* quadPtsRef; ///< Coordinates of quad pts in ref cell
  double* quadWts; ///< Weights of quadrature points
  double* quadPts; ///< Coordinates of quad pts in cell

  double* basis; ///< Basis fns at quadrature points
  double* basisDerivRef; ///< Derivatices of basis fns at quad pts (cell)
  double* basisDeriv; ///< Derivatices of basis fns at quad pts (global)
  double* jacobian; ///< Jacobian at quadrature points
  double* jacobianDet; ///< Determinant of quadrature points
  double* jacobianInv; ///< Inverse of Jacobian at quadruature points

};

#endif // pylith_feassemble_quadraturedata_hh

// End of file
