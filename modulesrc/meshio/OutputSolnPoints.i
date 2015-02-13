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

/**
 * @file modulesrc/meshio/OutputSolnPoints.i
 *
 * @brief Python interface to C++ OutputSolnPoints object.
 */

namespace pylith {
  namespace meshio {

    class pylith::meshio::OutputSolnPoints : public OutputManager
    { // OutputSolnPoints
      
      // PUBLIC METHODS ///////////////////////////////////////////////////
    public :
      
      /// Constructor
      OutputSolnPoints(void);
      
      /// Destructor
      ~OutputSolnPoints(void);
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
      
      /** Get mesh associated with points.
       *
       * @returns Mesh associated with points.
       */
      const pylith::topology::Mesh& pointsMesh(void);

      /** Setup interpolator.
       *
       * @param mesh Domain mesh.
       * @param points Array of dimensioned coordinates for points [numPoints*spaceDim].
       * @param numPoints Number of points.
       * @param spaceDim Spatial dimension for coordinates.
       */
      %apply(PylithScalar* IN_ARRAY2, int DIM1, int DIM2) {
	(const PylithScalar* points,
	 const int numPoints,
	 const int spaceDim)
      };
      void setupInterpolator(pylith::topology::Mesh* mesh,
			     const PylithScalar* points,
			     const int numPoints,
			     const int spaceDim,
			     const spatialdata::units::Nondimensional& normalizer);
      %clear(const PylithScalar* points, const int numPoints, const int spaceDim);
  
      /** Prepare for output.
       *
       * @param mesh Finite-element mesh object.
       * @param numTimeSteps Expected number of time steps.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void open(const pylith::topology::Mesh& mesh,
		const int numTimeSteps,
		const char* label =0,
		const int labelId =0);

      /** Setup file for writing fields at time step.
       *
       * @param t Time of time step.
       * @param mesh Finite-element mesh object.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void openTimeStep(const PylithScalar t,
			const pylith::topology::Mesh& mesh,
			const char* label =0,
			const int labelId =0);

      /** Append finite-element vertex field to file.
       *
       * @param t Time associated with field.
       * @param field Vertex field.
       * @param mesh Mesh for output.
       */
      void appendVertexField(const PylithScalar t,
			     pylith::topology::Field& field,
			     const pylith::topology::Mesh& mesh);

      /** Append finite-element cell field to file.
       *
       * @param t Time associated with field.
       * @param field Cell field.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void appendCellField(const PylithScalar t,
			   pylith::topology::Field& field,
			   const char* label =0,
			   const int labelId =0);

      /** Write dataset with names of points to file.
       *
       * @param names Array with name for each point, e.g., station name.
       * @param nunNames Number of names in array.
       */
      %apply(const char* const* string_list, const int list_len){
	(const char* const* names, const int numNames)
	  };
      void writePointNames(const char* const* names,
			   const int numNames);
      %clear(const char* const* names, const int numNames);

    }; // OutputSolnPoints

  } // meshio
} // pylith


// End of file 
