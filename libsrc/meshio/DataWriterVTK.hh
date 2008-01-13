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
   */
  void openTimeStep(const double t,
		    const ALE::Obj<ALE::Mesh>& mesh,
		    const spatialdata::geocoords::CoordSys* csMesh);

  /// Cleanup after writing data for a time step.
  void closeTimeStep(void);

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over vertices.
   * @param mesh Finite-element mesh
   */
  void writeVertexField(const double t,
			const char* name,
			const ALE::Obj<real_section_type>& field,
			const ALE::Obj<ALE::Mesh>& mesh);

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over cells.
   * @param mesh PETSc mesh object.
   */
  void writeCellField(const double t,
		      const char* name,
		      const ALE::Obj<real_section_type>& field,
		      const ALE::Obj<ALE::Mesh>& mesh);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _filename; ///< Name of VTK file.
  std::string _timeFormat; ///< C style time format for time stamp.

  PetscViewer _viewer; ///< Output file

}; // DataWriterVTK

#include "DataWriterVTK.icc" // inline methods

#endif // pylith_meshio_datawritervtk_hh


// End of file 
