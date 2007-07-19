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

/**
 * @file pylith/feassemble/CellGeometry.hh
 *
 * @brief C++ abstract base class for cell geometry calculations.
 */

#if !defined(pylith_feassemble_cellgeometry_hh)
#define pylith_feassemble_cellgeometry_hh

#include "pylith/utils/arrayfwd.hh" // USES double_array

namespace pylith {
  namespace feassemble {
    class CellGeometry;

    class TestCellGeometry; // unit testing
  } // feassemble
} // pylith

class pylith::feassemble::CellGeometry
{ // CellGeometry
  friend class TestCellGeometry; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /** Default constructor.
   *
   * @param cellDim Dimension of cell.
   * @param spaceDim Dimension of coordinate space.
   * @param numCorners Number of corners in cell.
   */
  CellGeometry(const int cellDim,
	       const int spaceDim,
	       const int numCorners);

  /// Default destructor.
  virtual
  ~CellGeometry(void);

  /** Create a copy of geometry.
   *
   * @returns Copy of geometry.
   */
  virtual
  CellGeometry* clone(void) const = 0;

  /** Get dimension of cell.
   *
   * @returns Spatial dimension of cell.
   */
  int cellDim(void) const;

  /** Get dimension of coordinate space.
   *
   * @returns Dimension of coordinate space.
   */
  int spaceDim(void) const;

  /** Get number of vertices in cell.
   *
   * @returns Number of vertices in cell.
   */
  int numCorners(void) const;

  /** Get cell geometry for lower dimension cell.
   *
   * @returns Pointer to cell geometry object corresponding to next
   * lower dimension, NULL if there is no lower dimension object.
   */
  virtual
  CellGeometry* geometryLowerDim(void) const = 0;

  /** Compute Jacobian at location in cell.
   *
   * @param jacobian Jacobian at location.
   * @param det Determinant of Jacobian at location.
   * @param vertices Coordinates of vertices of cell.
   * @param location Location in reference cell at which to compute Jacobian.
   */
  virtual
  void jacobian(double_array* jacobian,
		double* det,
		const double_array& vertices,
		const double_array& location) const = 0;

  /** Compute orientation of cell at location.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = spaceDim*spaceDim
   * index = iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to horizontal direction that is 
   *   not collinear with cell normal (usually "up" direction).
   */
  void orientation(double_array* orientation,
		   const double_array& jacobian,
		   const double jacobianDet,
		   const double_array& upDir) const;

// PROTECTED ////////////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param g Geometry to copy.
   */
  CellGeometry(const CellGeometry& g);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented
  const CellGeometry& operator=(const CellGeometry& );

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

  /// Function type for orientation methods.
  typedef void (*orient_fn_type)(double_array*, 
				 const double_array&,
				 const double,
				 const double_array&);  

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Compute weighted orientation of fault for cohesive cell between
   * 1-D elements. Orientation is either at vertices or quadrature
   * points, depending on whether the arguments have been evaluated at
   * the vertices or quadrature points.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = spaceDim*spaceDim
   * index = iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  static
  void _orient1D(double_array* orientation,
		 const double_array& jacobian,
		 const double jacobianDet,
		 const double_array& upDir);
		
  /** Compute weighted orientation of fault for cohesive cell between
   * 2-D elements. Orientation is either at vertices or quadrature
   * points, depending on whether the arguments have been evaluated at
   * the vertices or quadrature points.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = spaceDim*spaceDim
   * index = iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  static 
  void _orient2D(double_array* orientation,
		 const double_array& jacobian,
		 const double jacobianDet,
		 const double_array& upDir);
		
  /** Compute weighted orientation of fault for cohesive cell between
   * 3-D elements. Orientation is either at vertices or quadrature
   * points, depending on whether the arguments have been evaluated at
   * the vertices or quadrature points.
   *
   * The orientation is returned as an array of direction cosines.
   *
   * size = spaceDim*spaceDim
   * index = iDir*spaceDim + iComponent
   *
   * @param orientation Array of direction cosines.
   * @param jacobian Jacobian matrix at point.
   * @param jacobianDet Determinant of Jacobian matrix at point.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction).
   */
  static
  void _orient3D(double_array* orientation,
		 const double_array& jacobian,
		 const double jacobianDet,
		 const double_array& upDir);
		
// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  orient_fn_type _orientFn; ///< Function for computing orientation of cell  
  int _cellDim; ///< Dimension of cell.
  int _spaceDim; ///< Dimension of coordinate space.
  int _numCorners; ///< Number of corners in cell.

}; // CellGeometry

#include "CellGeometry.icc" // inline methods

#endif // pylith_feassemble_cellgeometry_hh


// End of file
