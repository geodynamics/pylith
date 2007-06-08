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

#if !defined(pylith_meshio_geomdatatri2d_hh)
#define pylith_meshio_geomdatatri2d_hh

#include "CellGeomData.hh"

namespace pylith {
  namespace feassemble {
     class GeomDataTri2D;
  } // feassemble
} // pylith

class pylith::feassemble::GeomDataTri2D : public CellGeomData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public: 

  /// Constructor
  GeomDataTri2D(void);

  /// Destructor
  ~GeomDataTri2D(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

  static const int _cellDim; ///< Number of dimensions associated with cell
  static const int _spaceDim; ///< Number of dimensions in vertex coordinates
  static const int _numCorners; ///< Number of vertices in cell

  static const int _numLocs; ///< Number of locations for computing Jacobian

  static const double _vertices[]; ///< Coordinates of cell's vertices
  static const double _locations[]; ///< Locations to compute Jacobian
  static const double _jacobian[]; ///< Jacobian at locations
  static const double _jacobianDet[]; ///< Determinant of Jacobian at locations

};

#endif // pylith_meshio_geomdatatri2d_hh

// End of file
