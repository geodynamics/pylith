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

  int numDOF; ///< Number of degrees of freedom at each point.
  int boundaryCellDim; ///< Dimension of surface cells.
  int numBoundaryCells; ///< Number of cells on Neumann boundary.
  int numBoundaryVertices; ///< Number of vertices on Neumann boundary.
  int numVertices; ///< Number of vertices in the mesh.
  int spaceDim; ///< Spatial dimension for the problem.

  int id; ///< Boundary condition identifier
  char* label; ///< Label for boundary condition group

  int* numCorners; ///< Array defining the number of vertices for each
                   /// boundary cell.
  int* cells; ///< Array of vertices defining each boundary cell.
  double* tractionVals; ///< Traction values at specified points.

  char* meshFilename; ///< Filename for input mesh.
  char* dbFilename; ///< Filename of simple spatial database.
};

#endif // pylith_bc_neumanndata_hh

// End of file
