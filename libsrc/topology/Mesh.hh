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
#include "topologyfwd.hh" // forward declarations
#include "spatialdata/geocoords/geocoordsfwd.hh" // forward declarations

#include <petscmesh.hh> // HASA ALE::IMesh

// Mesh -----------------------------------------------------------------
class pylith::topology::Mesh
{ // Mesh
  friend class TestMesh; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  // Typedefs for basic types associated with Sieve mesh.
  // All other PyLith mesh and submesh objects should define:
  //   (1) SieveMesh - Sieve mesh
  //   (2) RealSection - Section of doubles
  //   (3) IntSection - Section of ints
  // because these are used in templated code.
  // 
  // All other mesh objects for the domain should also define
  //   (1) SieveSubMesh - SubMesh object
  typedef ALE::IMesh<> SieveMesh;
  typedef SieveMesh::real_section_type RealSection;
  typedef SieveMesh::int_section_type IntSection;
  typedef ALE::IMesh<ALE::LabelSifter<int, SieveMesh::point_type> > SieveSubMesh;

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
