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
 * @file modulesrc/meshio/DataWriterHDF5Ext.i
 *
 * @brief Python interface to C++ DataWriterHDF5Ext object.
 */

namespace pylith {
  namespace meshio {

    template<typename mesh_type, typenam field_type>
    class pylith::meshio::DataWriterHDF5Ext :
      public DataWriter<mesh_type, field_type>
    { // DataWriterHDF5Ext  
      
      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      DataWriterHDF5Ext(void);
      
      /// Destructor
      ~DataWriterHDF5Ext(void);
      
      /** Make copy of this object.
       *
       * @returns Copy of this.
       */
      DataWriter<mesh_type, field_type>* clone(void) const;
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set filename for HDF5Ext file.
       *
       * @param filename Name of HDF5Ext file.
       */
      void filename(const char* filename);
      
      /** Open output file.
       *
       * @param mesh Finite-element mesh. 
       * @param numTimeSteps Expected number of time steps for fields.
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

      /** Write field over vertices to file.
       *
       * @param t Time associated with field.
       * @param field Field over vertices.
       * @param mesh Mesh for output.
       */
      void writeVertexField(const double t,
			    field_type& field,
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
			  field_type& field,
			  const char* label =0,
			  const int labelId =0);
      
    }; // DataWriterHDF5Ext

  } // meshio
} // pylith


// End of file 
