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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/meshio/DataWriterHDF5.i
 *
 * @brief Python interface to C++ DataWriterHDF5 object.
 */

namespace pylith {
  namespace meshio {

    template<typename mesh_type, typenam field_type>
    class pylith::meshio::DataWriterHDF5 :
      public DataWriter<mesh_type, field_type>
    { // DataWriterHDF5  
      
      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      DataWriterHDF5(void);
      
      /// Destructor
      ~DataWriterHDF5(void);
      
      /** Make copy of this object.
       *
       * @returns Copy of this.
       */
      DataWriter<mesh_type, field_type>* clone(void) const;
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set filename for HDF5 file.
       *
       * @param filename Name of HDF5 file.
       */
      void filename(const char* filename);
      
      /** Prepare file for data at a new time step.
       *
       * @param t Time stamp for new data
       * @param mesh Finite-element mesh.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void openTimeStep(const double t,
			const mesh_type& mesh,
			const char* label =0,
			const int labelId =0);
      
      /// Cleanup after writing data for a time step.
      void closeTimeStep(void);
      
      /** Write field over vertices to file.
       *
       * @param t Time associated with field.
       * @param field Field over vertices.
       * @param mesh Mesh for output.
       */
      void writeVertexField(const double t,
			    const field_type& field,
			    const mesh_type& mesh);
      
      /** Write field over cells to file.
       *
       * @param t Time associated with field.
       * @param field Field over cells.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void writeCellField(const double t,
			  const field_type& field,
			  const char* label =0,
			  const int labelId =0);
      
    }; // DataWriterHDF5

  } // meshio
} // pylith


// End of file 
