// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/** @file modulesrc/feassemble/CellGeometry.i
 *
 * @brief Python interface to C++ CellGeometry object.
 */

namespace pylith {
  namespace feassemble {

    class CellGeometry
    { // CellGeometry

      // PUBLIC ENUMS ///////////////////////////////////////////////////
    public :
      
      enum ShapeEnum { 
	LINE=2, // 1-D line cell (2 points)
	TRIANGLE=4, // 2-D triangular cell (3 edges)
	QUADRILATERAL=5, // 2-D quadrilateral cell (4 edges)
	TETRAHEDRON=16, // 3-D tetrahedral cell (4 faces)
	HEXAHEDRON=17 // 3-D hexahedral cell (6 faces)
      };

      // PUBLIC METHODS /////////////////////////////////////////////////
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
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
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
      const pylith::scalar_array& vertices(void) const;
      
      /** Get cell geometry for lower dimension cell.
       *
       * @returns Pointer to cell geometry object corresponding to next
       * lower dimension, NULL if there is no lower dimension object.
       */
      virtual
      CellGeometry* geometryLowerDim(void) const = 0;
      
      /** Transform coordinates in reference cell to global coordinates.
       *
       * @param ptsGlobal Array of points in global coordinate system.
       * @param ptsRef Array of points in reference cell.
       * @param vertices Array of cell vertices in global coordinates.
       * @param dim Dimension of global coordinate system.
       * @param npts Number of points to transform.
       */
      virtual
      void ptsRefToGlobal(PylithScalar* ptsGlobal,
			  const PylithScalar* ptsRef,
			  const PylithScalar* vertices,
			  const int dim,
			  const int npts) const = 0;
      
      /** Compute Jacobian at location in cell.
       *
       * @param jacobian Jacobian at location.
       * @param det Determinant of Jacobian at location.
       * @param vertices Coordinates of vertices of cell.
       * @param numVertices Number of vertices in cell.
       * @param spaceDim Spatial dimension of coordinates.
       * @param location Location in reference cell at which to compute Jacobian.
       * @param cellDim Dimension of reference cell.
       */
      virtual
      void jacobian(pylith::scalar_array* jacobian,
		    PylithScalar* det,
		    const PylithScalar* vertices,
		    const int numVertices,
		    const int spaceDim,
		    const PylithScalar* location,
		    const int cellDim) const = 0;

      /** Compute Jacobian at location in cell.
       *
       * @param jacobian Jacobian at location.
       * @param det Determinant of Jacobian at location.
       * @param vertices Coordinates of vertices of cell.
       * @param ptsRef Points in reference cell at which to compute Jacobian.
       * @param dim Dimension of coordinate system.
       * @param npts Number of points to transform.
       */
      virtual
      void jacobian(PylithScalar* jacobian,
		    PylithScalar* det,
		    const PylithScalar* vertices,
		    const PylithScalar* ptsRef,
		    const int dim,
		    const int npts) const = 0;
      
      /** Compute orientation of cell at location.
       *
       * The orientation is returned as an array of direction cosines
       * weighted by the determinant of the Jacobian.
       *
       * size = spaceDim*spaceDim
       * index = iDir*spaceDim + iComponent
       *
       * @param orientation Array of direction cosines (weighted by Jacobian).
       * @param jacobian Jacobian matrix at point.
       * @param jacobianDet Determinant of Jacobian matrix at point.
       * @param upDir Direction perpendicular to horizontal direction that is 
       *   not collinear with cell normal (usually "up" direction).
       */
      void orientation(pylith::scalar_array* orientation,
		       const pylith::scalar_array& jacobian,
		       const PylithScalar jacobianDet,
		       const pylith::scalar_array& upDir) const;

    }; // CellGeometry

  } // feassemble
} // pylith


// End of file
