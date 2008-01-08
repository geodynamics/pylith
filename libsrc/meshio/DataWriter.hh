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
 * @file pylith/meshio/DataWriter.hh
 *
 * @brief Abstract base class for writing finite-element data to file.
 */

#if !defined(pylith_meshio_datawriter_hh)
#define pylith_meshio_datawriter_hh

#include "pylith/utils/petscfwd.h" // USES PetscVec
#include "pylith/utils/sievetypes.hh" // USES ALE::Mesh, real_section_type

namespace pylith {
  namespace meshio {
    class DataWriter;
  } // meshio
} // pylith

namespace spatialdata {
  namespace geocoords {
    class CoordSys; // HOLDSA CoordSys
  } // geocoords
} // spatialdata

class pylith::meshio::DataWriter
{ // DataWriter

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  DataWriter(void);

  /// Destructor
  virtual
  ~DataWriter(void);

  /** Set coordinate system for output.
   *
   * @param cs Coordinate system
   */
  void coordsys(const spatialdata::geocoords::CoordSys* cs);

  /** Prepare for writing files.
   *
   * @param mesh PETSc mesh object 
   * @param csMesh Coordinate system of mesh geometry
   */
  virtual
  void open(const ALE::Obj<ALE::Mesh>& mesh,
	    const spatialdata::geocoords::CoordSys* csMesh);

  /// Close output files.
  virtual
  void close(void);

  /** Prepare file for data at a new time step.
   *
   * @param t Time stamp for new data
   * @param mesh PETSc mesh object
   * @param csMesh Coordinate system of mesh geometry
   */
  virtual
  void openTimeStep(const double t,
		    const ALE::Obj<ALE::Mesh>& mesh,
		    const spatialdata::geocoords::CoordSys* csMesh);

  /// Cleanup after writing data for a time step.
  virtual
  void closeTimeStep(void);

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param vec PETSc Vec field over vertices.
   * @param name Name of field.
   */
  virtual
  void writeVertexField(const double t,
			const PetscVec* vec,
			const char* name) = 0;

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param vec PETSc Vec field over cells.
   * @param name Name of field.
   */
  virtual
  void writeCellField(const double t,
		      const PetscVec* vec,
		      const char* name) = 0;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
public :

  spatialdata::geocoords::CoordSys* _cs; ///< Coordinate system for output

}; // DataWriter

#endif // pylith_meshio_datawriter_hh


// End of file 
