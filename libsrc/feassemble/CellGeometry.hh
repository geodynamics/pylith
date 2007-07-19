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

#include "pylith/utils/array.hh" // HASA double_array

namespace pylith {
  namespace feassemble {
    class CellGeometry;

    class TestCellGeometry; // unit testing
  } // feassemble
} // pylith

class pylith::feassemble::CellGeometry
{ // CellGeometry
  friend class TestCellGeometry; // unit testing

// PUBLIC ENUMS /////////////////////////////////////////////////////////
public :
  
  enum ShapeEnum { 
    POINT=0, // 0-D point cell
    LINE=2, // 1-D line cell (2 points)
    TRIANGLE=4, // 2-D triangular cell (3 edges)
    QUADRILATERAL=5, // 2-D quadrilateral cell (4 edges)
    TETRAHEDRON=16, // 3-D tetrahedral cell (4 faces)
    HEXAHEDRON=17 // 3-D hexahedral cell (6 faces)
  };

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :
  
  /** Default constructor.
   *
   * @param shape String identifying shape of cell
   * @param spaceDim Dimension of coordinate space.
   */
  CellGeometry(const ShapeEnum shape,
	       const int spaceDim);

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

  /** Get coordinates of vertices in reference cell (dual basis).
   *
   * @returns Array of coordinates of vertices in reference cell
   */
  const double_array& vertices(void) const;

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

  /** Set coordinates of vertices in reference cell.
   *
   * @param vertices Array of coordinates of vertices [#vertices*cellDim]
   * @param numVertices Number of vertices.
   * @param dim Dimension of cell
   */
  void _setVertices(const double* vertices,
		    const int numVertices,
		    const int dim);

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

  /** Compute orientation of 0-D cell. Orientation is either at
   * vertices or quadrature points, depending on whether the arguments
   * have been evaluated at the vertices or quadrature points.
   *
   * The orientation of a 0-D cell is always [1.0].
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
  void _orient0D(double_array* orientation,
		 const double_array& jacobian,
		 const double jacobianDet,
		 const double_array& upDir);
		
  /** Compute orientation of 1-D cell. Orientation is either at
   * vertices or quadrature points, depending on whether the arguments
   * have been evaluated at the vertices or quadrature points.
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
		
  /** Compute orientation of 2-D cell. Orientation is either at
   * vertices or quadrature points, depending on whether the arguments
   * have been evaluated at the vertices or quadrature points.
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
		
// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /** Array of coordinates of vertices in reference cell (dual basis).
   *
   * Reference coordinates: (p,q,r)
   *
   * v0p, v0q, v0r
   * v1p, v1q, v1r
   *
   * size = numBasis * cellDim
   * index = iBasis*cellDim + iDim
   */
  double_array _vertices;

  orient_fn_type _orientFn; ///< Function for computing orientation of cell
  int _spaceDim; ///< Dimension of coordinate space.
  ShapeEnum _shape; ///< Shape of cell

}; // CellGeometry

#include "CellGeometry.icc" // inline methods

#endif // pylith_feassemble_cellgeometry_hh


// End of file
