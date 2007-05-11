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

  /** Write solution topology to file.
   *
   * @param mesh PETSc mesh object.
   * @param 
   */
  virtual
  void writeTopology(const ALE::Obj<ALE::Mesh>& mesh,
		     const spatialdata::geocoords::CoordSys* csMesh) = 0;

  /** Write solution field to file.
   *
   * @param t Time associated with field.
   * @param field PETSc field.
   * @param name Name of field.
   * @param mesh PETSc mesh object.
   */
  virtual
  void writeField(const double t,
		  const ALE::Obj<real_section_type>& field,
		  const char* name,
		  const ALE::Obj<ALE::Mesh>& mesh) = 0;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
public :

  spatialdata::geocoords::CoordSys* _cs; ///< Coordinate system for output

}; // SolutionIO

#endif // pylith_meshio_solutionio_hh

// End of file 
