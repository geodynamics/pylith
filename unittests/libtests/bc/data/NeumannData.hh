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

#if !defined(pylith_bc_neumanndata_hh)
#define pylith_bc_neumanndata_hh

namespace pylith {
  namespace bc {
     class NeumannData;
  } // pylith
} // bc

class pylith::bc::NeumannData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  NeumannData(void);

  /// Destructor
  ~NeumannData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  char* meshFilename; ///< Filename for input mesh.

  /// @name Quadrature information
  //@{
  int spaceDim; ///< Spatial dimension for the problem.
  int cellDim; ///< Dimension of surface cells.
  int numBasis; ///< Number of basis functions for surface cells.
  int numQuadPts; ///< Number of quadrature points per boundary cell.
  double* quadPts; ///< Coordinates of quadrature points in ref cell.
  double* quadWts; ///< Weights of quadrature points.
  double* basis; ///< Cell basis functions at quad points.
  double* basisDeriv; ///< Derivatives of cell basis functions at quad points.
  double* verticesRef; ///< Coordinates of vertices in ref cell (dual basis).
  //@}

  /// @name BC information
  //@{
  int id; ///< Boundary condition identifier
  char* label; ///< Label for boundary condition group
  char* dbFilename; ///< Filename of simple spatial database.
  //@}

  /// @name Calculated values.
  //@{
  int numBoundaryCells; ///< Expected number of cells on Neumann boundary.
  int numVertices; ///< Expected number of vertices in the mesh.
  int* numCorners; ///< Expected number of vertices for each boundary cell.
  int* cells; ///< Expected array of vertices defining each boundary cell.
  double* tractionCell; ///< Expected traction values at quadrature points.
  double* valsResidual; ///< Expected residual at each vertex.
  //@}


};

#endif // pylith_bc_neumanndata_hh

// End of file
