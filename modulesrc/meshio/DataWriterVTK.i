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
 * @file modulesrc/meshio/DataWriterVTK.i
 *
 * @brief Python interface to C++ DataWriterVTK object.
 */

namespace pylith {
  namespace meshio {

    class pylith::meshio::DataWriterVTK : public DataWriter
    { // DataWriterVTK  
      
      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      DataWriterVTK(void);
      
      /// Destructor
      ~DataWriterVTK(void);
      
      /** Make copy of this object.
       *
       * @returns Copy of this.
       */
      DataWriter* clone(void) const;
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);
  
      /** Set filename for VTK file.
       *
       * @param filename Name of VTK file.
       */
      void filename(const char* filename);
      
      /** Set time format for time stamp in name of VTK file.
       *
       * @param format C style time format for filename.
       */
      void timeFormat(const char* format);
      
      /** Set value used to normalize time stamp in name of VTK file.
       *
       * Time stamp is divided by this value (time in seconds).
       *
       * @param value Value (time in seconds) used to normalize time stamp in
       * filename.
       */
      void timeConstant(const PylithScalar value);
      
      /** Set precision of floating point values in output.
       *	
       * @param value Precision for floating point values.
       */
      void precision(const int value);

      /** Prepare for writing files.
       *
       * @param mesh Finite-element mesh. 
       * @param numTimeSteps Expected number of time steps for fields.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void open(const pylith::topology::Mesh& mesh,
		const int numTimeSteps,
		const char* label =0,
		const int labelId =0);
      
      /// Close output files.
      void close(void);

      /** Prepare file for data at a new time step.
       *
       * @param t Time stamp for new data
       * @param mesh Finite-element mesh.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void openTimeStep(const PylithScalar t,
			const pylith::topology::Mesh& mesh,
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
      void writeVertexField(const PylithScalar t,
			    pylith::topology::Field& field,
			    const pylith::topology::Mesh& mesh);
      
      /** Write field over cells to file.
       *
       * @param t Time associated with field.
       * @param field Field over cells.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void writeCellField(const PylithScalar t,
			  pylith::topology::Field& field,
			  const char* label =0,
			  const int labelId =0);
      
    }; // DataWriterVTK

  } // meshio
} // pylith


// End of file 
