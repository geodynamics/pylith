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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/meshio/OutputManager.i
 *
 * @brief Python interface to C++ OutputManager object.
 */

namespace pylith {
  namespace meshio {

    template<typename mesh_type, typename field_type>
    class pylith::meshio::OutputManager
    { // OutputManager

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      OutputManager(void);
      
      /// Destructor
      virtual
      ~OutputManager(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set coordinate system in output. The vertex fields in the output
       * are not affected by any change in coordinates.
       *
       * @param cs Coordinate system in output.
       */
      void coordsys(const spatialdata::geocoords::CoordSys* cs);
      
      /** Set writer to write data to file.
       *
       * @param datawriter Writer for data.
       */
      void writer(DataWriter<mesh_type, field_type>* const datawriter);
      
      /** Set filter for vertex data.
       *
       * @param filter Filter to apply to vertex data before writing.
       */
      void vertexFilter(VertexFilter<field_type>* const filter);
      
      /** Set filter for cell data.
       *
       * @param filter Filter to apply to cell data before writing.
       */
      void cellFilter(CellFilter<mesh_type, field_type>* const filter);
      
      /** Get fields used in output.
       *
       * @returns Fields associated with output.
       */
      const pylith::topology::Fields<field_type>* fields(void) const;

      /** Prepare for output.
       *
       * @param mesh Finite-element mesh object.
       * @param numTimeSteps Expected number of time steps.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void open(const mesh_type& mesh,
		const int numTimeSteps,
		const char* label =0,
		const int labelId =0);
      
      /// Close output files.
      void close(void);
      
      /** Setup file for writing fields at time step.
       *
       * @param t Time of time step.
       * @param mesh Finite-element mesh object.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void openTimeStep(const double t,
			const mesh_type& mesh,
			const char* label =0,
			const int labelId =0);
      
      /// End writing fields at time step.
      void closeTimeStep(void);
      
      /** Append finite-element vertex field to file.
       *
       * @param t Time associated with field.
       * @param field Vertex field.
       * @param mesh Mesh for output.
       */
      void appendVertexField(const double t,
			     field_type& field,
			     const mesh_type& mesh);
      
      /** Append finite-element cell field to file.
       *
       * @param t Time associated with field.
       * @param field Cell field.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void appendCellField(const double t,
			   field_type& field,
			   const char* label =0,
			   const int labelId =0);

    }; // OutputManager

  } // meshio
} // pylith


// End of file 
