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

#if !defined(pylith_bc_neumanndatahex8_hh)
#define pylith_bc_neumanndatahex8_hh

#include "NeumannData.hh"

namespace pylith {
  namespace bc {
     class NeumannDataHex8;
  } // pylith
} // bc

class pylith::bc::NeumannDataHex8 : public NeumannData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  NeumannDataHex8(void);

  /// Destructor
  ~NeumannDataHex8(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  // Quadrature information
  static const char* _meshFilename; ///< Filename of input mesh.
  static const int _spaceDim; ///< Dimension of mesh.
  static const int _cellDim; ///< Dimension of surface cells.
  static const int _numBasis; ///< Number of basis functions for surface cells.
  static const int _numQuadPts; ///< Number of quadrature points per boundary cell.
  static const double _quadPts[]; ///< Coordinates of quadrature points in ref cell.
  static const double _quadWts[]; ///< Weights of quadrature points.
  static const double _basis[]; ///< Cell basis functions at quad points.
  static const double _basisDeriv[]; ///< Derivatives of cell basis functions at quad points.
  static const double _verticesRef[]; ///< Coordinates of vertices in ref cell (dual basis).

  // BC information
  static const int _id; ///< Boundary condition identifier
  static const char* _label; /// Label for boundary condition group
  static const char* _dbFilename; ///< Filename of simple spatial database.

  // Calculated values.
  static const int _numBoundaryCells; ///< Expected number of cells on Neumann boundary.
  static const int _numVertices; ///< Expected number of vertices in the mesh.
  static const int _numCorners[]; ///< Expected number of vertices for each boundary cell.
  static const int _cells[]; ///< Expected array of vertices defining each boundary cell.
  static const double _tractionsCell[]; ///< Expected traction values at quadrature points.
  static const double _valsResidual[]; ///< Expected residual at each vertex.

};

#endif // pylith_bc_neumanndatahex8_hh

// End of file
