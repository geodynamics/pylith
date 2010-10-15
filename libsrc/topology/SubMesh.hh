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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/SubMesh.hh
 *
 * @brief C++ PyLith finite-element mesh.
 */

#if !defined(pylith_topology_submesh_hh)
#define pylith_topology_submesh_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations
#include "spatialdata/geocoords/geocoordsfwd.hh" // forward declarations

#include "Mesh.hh" // USES Mesh

// SubMesh -----------------------------------------------------------------
/** @brief C++ PyLith finite-element mesh.
 *
 * Extends Sieve mesh over subset of domain to include coordinate
 * system associated with domain. Also has functions to simply
 * creating submeshes from groups of vertices.
 */
class pylith::topology::SubMesh
{ // SubMesh
  friend class TestSubMesh; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public:

  // Sieve mesh for higher level domain (mesh, not submesh)
  typedef Mesh::SieveMesh DomainSieveMesh;

  // Typedefs for basic types associated with Sieve mesh.
  // SieveMesh, RealSection, and IntSection are used in templated code.
  typedef Mesh::SieveSubMesh SieveMesh;

  typedef Mesh::IntSection IntSection;
  typedef Mesh::RealSection RealSection;
  typedef Mesh::RealUniformSection  RealUniformSection;

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  SubMesh(void);

  /** Constructor with mesh and label for vertices marking boundary.
   *
   * @param mesh Finite-element mesh over domain.
   * @param label Label for vertices marking boundary.
   */
  SubMesh(const Mesh& mesh,
	  const char* label);

  /// Default destructor
  ~SubMesh(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Create Sieve mesh.
   *
   * @param mesh Finite-element mesh over domain.
   * @param label Label for vertices marking boundary.
   */
  void createSubMesh(const Mesh& mesh,
		     const char* label);

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

  /** Set coordinate system using mesh.
   *
   * @param mesh Finite-element mesh over domain.
   */
  void coordsys(const Mesh& mesh);

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
  void view(const char* label) const;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  ALE::Obj<SieveMesh> _mesh; ///< Sieve mesh.
  spatialdata::geocoords::CoordSys* _coordsys; ///< Coordinate system.
  bool _debug; ///< Debugging flag for mesh.
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  SubMesh(const SubMesh&); ///< Not implemented
  const SubMesh& operator=(const SubMesh&); ///< Not implemented

}; // SubMesh

#include "SubMesh.icc"

#endif // pylith_topology_submesh_hh


// End of file
