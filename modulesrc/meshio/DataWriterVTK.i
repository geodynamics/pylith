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
 * @file modulesrc/meshio/DataWriterVTK.i
 *
 * @brief Python interface to C++ DataWriterVTK object.
 */

namespace pylith {
  namespace meshio {

    template<typename mesh_type>
    class pylith::meshio::DataWriterVTK : public DataWriter<mesh_type>
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
      DataWriter<mesh_type>* clone(void) const;
      
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
      void timeConstant(const double value);
      
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
       */
      void writeVertexField(const double t,
			    const pylith::topology::Field<mesh_type>& field);
      
      /** Write field over cells to file.
       *
       * @param t Time associated with field.
       * @param field Field over cells.
       * @param label Name of label defining cells to include in output
       *   (=0 means use all cells in mesh).
       * @param labelId Value of label defining which cells to include.
       */
      void writeCellField(const double t,
			  const pylith::topology::Field<mesh_type>& field,
			  const char* label =0,
			  const int labelId =0);
      
    }; // DataWriterVTK

  } // meshio
} // pylith


// End of file 
