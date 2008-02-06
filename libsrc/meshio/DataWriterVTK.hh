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
 * @file pylith/meshio/DataWriterVTK.hh
 *
 * @brief Abstract base class for writing finite-element data to file.
 */

#if !defined(pylith_meshio_datawritervtk_hh)
#define pylith_meshio_datawritervtk_hh

#include "DataWriter.hh" // ISA DataWriter

namespace pylith {
  namespace meshio {
    class DataWriterVTK;
  } // meshio
} // pylith

class pylith::meshio::DataWriterVTK : public DataWriter
{ // DataWriterVTK

// PUBLIC METHODS ///////////////////////////////////////////////////////
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

  /** Prepare file for data at a new time step.
   *
   * @param t Time stamp for new data
   * @param mesh PETSc mesh object
   * @param csMesh Coordinate system of mesh geometry
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  void openTimeStep(const double t,
		    const ALE::Obj<ALE::Mesh>& mesh,
		    const spatialdata::geocoords::CoordSys* csMesh,
		    const char* label =0,
		    const int labelId =0);

  /// Cleanup after writing data for a time step.
  void closeTimeStep(void);

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over vertices.
   * @param fieldType Type of field.
   * @param mesh Finite-element mesh
   */
  void writeVertexField(const double t,
			const char* name,
			const ALE::Obj<real_section_type>& field,
			const FieldEnum fieldType,
			const ALE::Obj<ALE::Mesh>& mesh);

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over cells.
   * @param fieldType Type of field.
   * @param mesh PETSc mesh object.
   */
  void writeCellField(const double t,
		      const char* name,
		      const ALE::Obj<real_section_type>& field,
		      const FieldEnum fieldType,
		      const ALE::Obj<ALE::Mesh>& mesh);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Copy constructor.
   *
   * @param w Object to copy.
   */
  DataWriterVTK(const DataWriterVTK& w);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  const DataWriterVTK& operator=(const DataWriterVTK&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _filename; ///< Name of VTK file.
  std::string _timeFormat; ///< C style time format for time stamp.

  PetscViewer _viewer; ///< Output file

  bool _wroteVertexHeader; ///< True if wrote header for vertex data.
  bool _wroteCellHeader; ///< True if wrote header for cell data

}; // DataWriterVTK

#include "DataWriterVTK.icc" // inline methods

#endif // pylith_meshio_datawritervtk_hh


// End of file 
