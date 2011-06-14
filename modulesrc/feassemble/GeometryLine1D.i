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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/** @file modulesrc/feassemble/GeometryLine1D.i
 *
 * @brief Python interface to C++ GeometryLine1D object.
 */

namespace pylith {
  namespace feassemble {

    class GeometryLine1D : public CellGeometry
    { // GeometryLine1D

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :
  
      /// Default constructor.
      GeometryLine1D(void);
      
      /// Default destructor.
      ~GeometryLine1D(void);
      
      /** Create a copy of geometry.
       *
       * @returns Copy of geometry.
       */
      CellGeometry* clone(void) const;
      
      /** Get cell geometry for lower dimension cell.
       *
       * @returns Pointer to cell geometry object corresponding to next
       * lower dimension, NULL if there is no lower dimension object.
       */
      CellGeometry* geometryLowerDim(void) const;
      
      /** Transform coordinates in reference cell to global coordinates.
       *
       * @param ptsGlobal Array of points in global coordinate system.
       * @param ptsRef Array of points in reference cell.
       * @param vertices Array of cell vertices in global coordinates.
       * @param dim Dimension of global coordinate system.
       * @param npts Number of points to transform.
       */
      void ptsRefToGlobal(double* ptsGlobal,
			  const double* ptsRef,
			  const double* vertices,
			  const int dim,
			  const int npts =1) const;
      
      /** Compute Jacobian at location in cell.
       *
       * @param jacobian Jacobian at location.
       * @param det Determinant of Jacobian at location.
       * @param vertices Coordinates of vertices of cell.
       * @param location Location in reference cell at which to compute Jacobian.
       */
      void jacobian(pylith::double_array* jacobian,
		    double* det,
		    const pylith::double_array& vertices,
		    const pylith::double_array& location) const;
      
      /** Compute Jacobian at location in cell.
       *
       * @param jacobian Jacobian at location.
       * @param det Determinant of Jacobian at location.
       * @param vertices Coordinates of vertices of cell.
       * @param ptsRef Points in reference cell at which to compute Jacobian.
       * @param dim Dimension of coordinate system.
       * @param npts Number of points to transform.
       */
      void jacobian(double* jacobian,
		    double* det,
		    const double* vertices,
		    const double* ptsRef,
		    const int dim,
		    const int npts =1) const;

    }; // GeometryLine1D  
    
  } // feassemble
} // pylith
  

// End of file
