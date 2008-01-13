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
 * @file pylith/meshio/SolutionIO.hh
 *
 * @brief Abstract base class for writing finite-element solution to file.
 *
 * :TODO: Replace this implementation with one that extracts data from
 * Sieve and calls virtual functions to write data to file (similar to
 * interface for MeshIO), so that children of SolutionIO don't
 * replicate code. This will also make it easier to implement
 * exporters for various formats (HDF5, OpenDX, etc).
 */

#if !defined(pylith_meshio_solutionio_hh)
#define pylith_meshio_solutionio_hh

#include "pylith/utils/sievetypes.hh" // USES ALE::Obj, ALE::Mesh, real_section_type

namespace pylith {
  namespace meshio {
    class SolutionIO;
  } // meshio
} // pylith

namespace spatialdata {
  namespace geocoords {
    class CoordSys; // HOLDSA CoordSys
  } // geocoords
} // spatialdata

class pylith::meshio::SolutionIO
{ // SolutionIO

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  SolutionIO(void);

  /// Destructor
  virtual
  ~SolutionIO(void);

  /** Set coordinate system for output.
   *
   * @param cs Coordinate system
   */
  void coordsys(const spatialdata::geocoords::CoordSys* cs);

  /** Open output files.
   *
   * @param mesh PETSc mesh object
   */
  virtual
  void open(const ALE::Obj<ALE::Mesh>& mesh) = 0;

  /// Close output files.
  virtual
  void close(void) = 0;

  /** Prepare file for data at a new time step.
   *
   * @param t Time stamp for new data
   * @param mesh PETSc mesh object
   * @param csMesh Coordinate system of mesh geometry
   */
  virtual
  void openTimeStep(const double t,
		    const ALE::Obj<ALE::Mesh>& mesh,
		    const spatialdata::geocoords::CoordSys* csMesh) = 0;

  /// Cleanup after writing data for a time step.
  virtual
  void closeTimeStep(void) = 0;

  /** Write field over vertices to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over vertices.
   * @param mesh PETSc mesh object.
   */
  virtual
  void writeVertexField(const double t,
			const char* name,
			const ALE::Obj<real_section_type>& field,
			const ALE::Obj<ALE::Mesh>& mesh) = 0;

  /** Write field over cells to file.
   *
   * @param t Time associated with field.
   * @param name Name of field.
   * @param field PETSc field over cells.
   * @param mesh PETSc mesh object.
   */
  virtual
  void writeCellField(const double t,
		      const char* name,
		      const ALE::Obj<real_section_type>& field,
		      const ALE::Obj<ALE::Mesh>& mesh) = 0;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
public :

  spatialdata::geocoords::CoordSys* _cs; ///< Coordinate system for output

}; // SolutionIO

#endif // pylith_meshio_solutionio_hh

// End of file 
