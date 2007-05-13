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

#if !defined(pylith_topology_cellgeomdata_hh)
#define pylith_topology_cellgeomdata_hh

namespace pylith {
  namespace topology {
     class CellGeomData;
  } // pylith
} // topology

class pylith::topology::CellGeomData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /// Constructor
  CellGeomData(void);

  /// Destructor
  ~CellGeomData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

  int cellDim; ///< Number of dimensions associated with cell
  int spaceDim; ///< Number of dimensions in vertex coordinates
  int numCorners; ///< Number of vertices in cell

  int numLocs; ///< Number of locations

  double* vertices; ///< Coordinates of vertices of cell
  double* locations; ///< Locations where Jacobian is computed
  double* jacobian; ///< Jacobian at locations

};

#endif // pylith_topology_cellgeomdata_hh


// End of file
