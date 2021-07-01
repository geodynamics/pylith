// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/Distributor.hh
 *
 * @brief Object for managing distribution of mesh among processors.
 */

#if !defined(pylith_topology_distributor_hh)
#define pylith_topology_distributor_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/meshio/meshiofwd.hh" // USES DataWriter<Mesh>

// Distributor ----------------------------------------------------------
/// Distribute mesh among processors.
class pylith::topology::Distributor
{ // Distributor
  friend class TestDistributor; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  Distributor(void);

  /// Destructor
  ~Distributor(void);

  /** Distribute mesh among processors.
   *
   * @param newMesh Distributed mesh (result).
   * @param origMesh Mesh to distribute.
   * @param partitionerName Name of PETSc partitioner to use in distributing mesh.
   */
  static
  void distribute(topology::Mesh* const newMesh,
		  const topology::Mesh& origMesh,
		  const char* partitionerName);

  /** Write partitioning info for distributed mesh.
   *
   * @param writer Data writer for partition information.
   * @param mesh Distributed mesh.
   * @param cs Coordinate system for mesh.
   */
  static
  void write(meshio::DataWriter* const writer,
	     const topology::Mesh& mesh);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Distributor(const Distributor&); ///< Not implemented
  const Distributor& operator=(const Distributor&); ///< Not implemented

}; // Distributor

#endif // pylith_topology_distributor_hh


// End of file 
