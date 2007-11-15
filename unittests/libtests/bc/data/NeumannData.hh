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
  int numBasis; ///< Number of basis functions for surface cells.
  int numQuadPts; ///< Number of quadrature points per boundary cell.
  double* quadPts; ///< Coordinates of quadrature points in ref cell.
  double* quadWts; ///< Weights of quadrature points.
  double* basis; ///< Cell basis functions at quad points.
  double* basisDerivRef; ///< Derivatives of basis functions at quad points.
  //@}

  /// @name Parameter information
  //@{
  char* spatialDBFilename; ///< Filename of simple spatial database.
  int id; ///< Boundary condition identifier
  char* label; ///< Label for boundary condition group
  //@}

  /// @name Boundary mesh information
  //@{
  int spaceDim; ///< Number of dimensions of vertex coordinates
  int cellDim; ///< Dimension of surface cells.
  int numBoundaryVertices; ///< Number of boundary vertices in the mesh.
  int numBoundaryCells; ///< Number of cells on Neumann boundary.
  int numCorners; ///< Number of vertices for each boundary cell.
  double* cellVertices; ///< Vertex coordinates for boundary cells.
  //@}

  /// @name Calculated values.
  //@{
  double* tractionsCell; ///< Expected traction values at quadrature points.
  double* valsResidual; ///< Expected residual at each vertex.
  //@}


};

#endif // pylith_bc_neumanndata_hh

// End of file
