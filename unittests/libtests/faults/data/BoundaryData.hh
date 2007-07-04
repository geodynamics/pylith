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

#if !defined(pylith_faults_boundarydata_hh)
#define pylith_faults_boundarydata_hh

namespace pylith {
  namespace faults {
     class BoundaryData;
  } // pylith
} // faults

class pylith::faults::BoundaryData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  BoundaryData(void);

  /// Destructor
  ~BoundaryData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int numVertices; ///< Number of vertices
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int numCells; ///< Number of cells
  int cellDim; ///< Number of dimensions associated with cell

  double* vertices; ///< Pointer to coordinates of vertices
  int* numCorners; ///< Number of vertices in cell
  int* cells; ///< Pointer to indices of vertices in cells

  char* filename; ///< Filename for input mesh
};

#endif // pylith_faults_boundarydata_hh

// End of file
