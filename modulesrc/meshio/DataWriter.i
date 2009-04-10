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
 * @file modulesrc/meshio/DataWriter.i
 *
 * @brief Python interface to C++ DataWriter object.
 */

namespace pylith {
  namespace meshio {

    template<typename mesh_type>
    class pylith::meshio::DataWriter
    { // DataWriter

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Constructor
      DataWriter(void);

      /// Destructor
      virtual
      ~DataWriter(void);
      
      /** Make copy of this object.
       *
       * @returns Copy of this.
       */
      virtual
      DataWriter* clone(void) const = 0;
      
      /** Prepare for writing files.
       *
       * @param mesh Finite-element mesh. 
       * @param numTimeSteps Expected number of time steps for fields.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      virtual
      void open(const mesh_type& mesh,
		const int numTimeSteps,
		const char* label =0,
		const int labelId =0);
      
      /// Close output files.
      virtual
      void close(void);
      
      /** Prepare file for data at a new time step.
       *
       * @param t Time stamp for new data
       * @param mesh PETSc mesh object
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      virtual
      void openTimeStep(const double t,
			const mesh_type& mesh,
			const char* label =0,
			const int labelId =0);
      
      /// Cleanup after writing data for a time step.
      virtual
      void closeTimeStep(void);
      
      /** Write field over vertices to file.
       *
       * @param t Time associated with field.
       * @param field Field over vertices.
       */
      virtual
      void writeVertexField(const double t,
			    const pylith::topology::Field<mesh_type>& field) = 0;
      
      /** Write field over cells to file.
       *
       * @param t Time associated with field.
       * @param field Field over cells.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      virtual
      void writeCellField(const double t,
			  const pylith::topology::Field<mesh_type>& field,
			  const char* label =0,
			  const int labelId =0) = 0;

    }; // DataWriter

  } // meshio
} // pylith


// End of file 
