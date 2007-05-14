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

#if !defined(pylith_meshio_geomdatahex3d_hh)
#define pylith_meshio_geomdatahex3d_hh

#include "CellGeomData.hh"

namespace pylith {
  namespace topology {
     class GeomDataHex3D;
  } // topology
} // pylith

class pylith::topology::GeomDataHex3D : public CellGeomData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  GeomDataHex3D(void);

  /// Destructor
  ~GeomDataHex3D(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _cellDim; ///< Number of dimensions associated with cell
  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _numCorners; ///< Number of vertices in cell

  static const int _numLocs; ///< Number of locations for computing Jacobian

  static const double _vertices[]; ///< Coordinates of cell's vertices
  static const double _locations[]; ///< Locations to compute Jacobian
  static const double _jacobian[]; ///< Jacobian at locations

};

#endif // pylith_meshio_geomdatahex3d_hh

// End of file
