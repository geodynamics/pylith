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

/** @file modulesrc/feassemble/Quadrature.i
 *
 * @brief Python interface to C++ Quadrature object.
 */

namespace pylith {
  namespace feassemble {

    class QuadratureRefCell
    { // Quadrature

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      QuadratureRefCell(void);

      /// Destructor
      virtual
      ~QuadratureRefCell(void);

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set basis functions and their derivatives, and coordinates and
       *  weights of the quadrature points.
       *
       * @param basis Array of basis functions evaluated at quadrature pts
       *   N0Qp0, N1Qp0, ...
       *   N0Qp1, N1Qp1, ...
       *   ...
       *   size = numQuadPts * numBasis
       *   index = iQuadPt*numBasis + iBasis
       *
       * @param basisDerivRef Array of basis function derivaties evaluated at
       * quadrature pts, where derivatives are with respect to cell's
       * local coordinates.
       *   N0pQp0, N0qQp0, N0rQp0, N1pQp0, N1qQp0, N1rQp0, ... 
       *   N0pQp1, N0qQp1, N0rQp1, N1pQp1, N1qQp1, N1rQp1, ...
       *   ...
       *   size = numQuadPts * numBasis * cellDim
       *   index = iQuadPt*numBasis*cellDim + iBasis*cellDim + iDim
       *
       * @param quadPts Array of coordinates of quadrature points in 
       *   reference cell
       *   Qp0p, Qp0q, Qp0r
       *   Qp1p, Qp1q, Qp1r
       *   size = numQuadPts * numDims
       *   index = iQuadPt*numDims + iDim
       *
       * @param quadWts Array of weights of quadrature points
       *   WtQp0, WtQp1, ...
       *   index = iQuadPt
       *
       * @param cellDim Number of dimensions in reference cell
       * @param numBasis Number of basis functions for a cell
       * @param numQuadPts Number of quadrature points
       * @param spaceDim Number of dimensions in coordinates of cell vertices
       */
      %apply(PylithScalar* IN_ARRAY2, int DIM1, int DIM2) {
	(const PylithScalar* basis,
	 const int numQuadPts1,
	 const int numBasis1)
	  };
      %apply(PylithScalar* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {
	(const PylithScalar* basisDerivRef,
	 const int numQuadPts2,
	 const int numBasis2,
	 const int cellDim2)
	  };
      %apply(PylithScalar* IN_ARRAY2, int DIM1, int DIM2) {
	(const PylithScalar* quadPtsRef,
	 const int numQuadPts3,
	 const int cellDim3)
	  };
      %apply(PylithScalar* IN_ARRAY1, int DIM1) {
	(const PylithScalar* quadWts,
	 const int numQuadPts4)
	  };
      void initialize(const PylithScalar* basis,
		      const int numQuadPts1,
		      const int numBasis1,
		      const PylithScalar* basisDerivRef,
		      const int numQuadPts2,
		      const int numBasis2,
		      const int cellDim2,
		      const PylithScalar* quadPtsRef,
		      const int numQuadPts3,
		      const int cellDim3,
		      const PylithScalar* quadWts,
		      const int numQuadPts4,
		      const int spaceDim);
      %clear(const PylithScalar* basis, const int numQuadPts, const int numBasis);
      %clear(const PylithScalar* basisDerivRef, const int numQuadPts,
	     const int numBasis, const int spaceDim);
      %clear(const PylithScalar* quadPtsRef, const int numQuadPts,
	     const int cellDim);
      %clear(const PylithScalar* quadWts, const int numQuadPts);
      
      /** Set geometry associated with reference cell.
       *
       * @param geometry Geometry of reference cell.
       */
      void refGeometry(CellGeometry* const geometry);
      
      /** Get geometry associated with reference cell.
       *
       * @returns Geometry of reference cell.
       */
      const CellGeometry& refGeometry(void) const;
      
      /** Set minimum allowable determinant of Jacobian.
       *
       * @param tolerance Minimum allowable value for Jacobian
       */
      void minJacobian(const PylithScalar min);
      
      /** Get minimum allowable determinant of Jacobian.
       *
       * @returns Minimum allowable value for Jacobian
       */
      PylithScalar minJacobian(void) const;
      
      /** Get coordinates of quadrature points in reference cell.
       *
       * @returns Array of coordinates of quadrature points in reference cell.
       */
      const pylith::scalar_array& quadPtsRef(void) const;
      
      /** Get weights of quadrature points.
       *
       * @returns Weights of quadrature points
       */
      const pylith::scalar_array& quadWts(void) const;
      
      /** Get basis fns evaluated at quadrature points.
       *
       * @returns Array of basis fns evaluated at quadrature points
       */
      const pylith::scalar_array& basis(void) const;
      
      /** Get derivates of basis fns evaluated at quadrature points.
       *
       * @returns Array of derivates of basis fns evaluated at
       * quadrature points
       */
      const pylith::scalar_array& basisDerivRef(void) const;
      
      /** Get number of dimensions in reference cell.
       *
       * @returns Number of dimensions in reference cell
       */
      int cellDim(void) const;
      
      /** Get number of basis functions for cell.
       *
       * @returns Number of basis functions for cell
       */
      int numBasis(void) const;
      
      /** Get number of quadrature points.
       *
       * @returns Number of quadrature points
       */
      int numQuadPts(void) const;
      
      /** Get number of dimensions in coordinates of cell vertices.
       *
       * @returns Number of dimensions in coordinates of cell vertices
       */
      int spaceDim(void) const;
      
    }; // QuadratureRefCell

  } // feassemble
} // pylith


// End of file 
