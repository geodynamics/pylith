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
 * @file libsrc/topology/Mesh.hh
 *
 * @brief C++ PyLith finite-element mesh.
 *
 * Extends Sieve mesh to include coordinate system associated with
 * domain.
 */

#if !defined(pylith_topology_mesh_hh)
#define pylith_topology_mesh_hh

// Include directives ---------------------------------------------------
#define NEWPYLITHMESH 1 
#include "pylith/utils/sievetypes.hh"

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace topology {
    class Mesh;
    class TestMesh; // unit testing
  } // topology
} // pylith

namespace spatialdata {
  namespace geocoords {
    class CoordSys;
  } // geocoords
} // spatialdata

// Mesh -----------------------------------------------------------------
class pylith::topology::Mesh
{ // Mesh
  friend class TestMesh; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  Mesh(void);

  /** Constructor with dimension and communicator.
   *
   * @param dim Dimension associated with mesh cells.
   * @param comm MPI communicator for mesh.
   */
  Mesh(const int dim,
       const MPI_Comm& comm =PETSC_COMM_WORLD); 

  /// Default destructor
  ~Mesh(void);

  /** Create Sieve mesh.
   *
   * @param dim Dimension associated with mesh cells.
   */
  void createSieveMesh(const int dim=3); 

  /** Get Sieve mesh.
   *
   * @returns Sieve mesh.
   */
  const ALE::Obj<SieveMesh>& sieveMesh(void) const;

  /** Get Sieve mesh.
   *
   * @returns Sieve mesh.
   */
  ALE::Obj<SieveMesh>& sieveMesh(void);

  /** Set coordinate system.
   *
   * @param cs Coordinate system.
   */
  void coordsys(const spatialdata::geocoords::CoordSys* cs);

  /** Get coordinate system.
   *
   * @returns Coordinate system.
   */
  const spatialdata::geocoords::CoordSys* coordsys(void) const;

  /** Set debug flag.
   *
   * @param value Turn on debugging if true.
   */
   void debug(const bool value);

  /** Get debug flag.
   *
   * @param Get debugging flag.
   */
   bool debug(void) const;

  /** Get dimension of mesh.
   *
   * @returns Dimension of mesh.
   */
  int dimension(void) const;

  /** Set MPI communicator associated with mesh.
   *
   * @param value MPI communicator.
   */
  void comm(const MPI_Comm value);
    
  /** Get MPI communicator associated with mesh.
   *
   * @returns MPI communicator.
   */
  const MPI_Comm comm(void) const;
    
  /// Initialize the finite-element mesh.
  void initialize(void);

  /** Print mesh to stdout.
   *
   * @param label Label for mesh.
   */
  void view(const char* label);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  ALE::Obj<SieveMesh> _mesh; ///< Sieve mesh.
  spatialdata::geocoords::CoordSys* _coordsys; ///< Coordinate system.
  MPI_Comm _comm; ///< MPI communicator for mesh.
  bool _debug; ///< Debugging flag for mesh.
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Mesh(const Mesh&); ///< Not implemented
  const Mesh& operator=(const Mesh&); ///< Not implemented

}; // Mesh

#include "Mesh.icc"

#endif // pylith_topology_mesh_hh


// End of file
