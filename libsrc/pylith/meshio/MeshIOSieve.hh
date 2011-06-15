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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/MeshIOSieve.hh
 *
 * @brief C++ input/output manager for Sieve mesh files.
 */

#if !defined(pylith_meshio_meshiosieve_hh)
#define pylith_meshio_meshiosieve_hh

// Include directives ---------------------------------------------------
#include "MeshIO.hh" // ISA MeshIO

// MeshIOSieve ----------------------------------------------------------
/// C++ input/output manager for Sieve mesh files.
class pylith::meshio::MeshIOSieve : public MeshIO
{ // MeshIOSieve
  friend class TestMeshIOSieve; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  MeshIOSieve(void);

  /// Destructor
  ~MeshIOSieve(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set filename for Sieve mesh file.
   *
   * @param filename Name of file
   */
  void filename(const char* name);

  /** Get filename of Sieve mesh file.
   *
   * @returns Name of file
   */
  const char* filename(void) const;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /// Write mesh
  void _write(void) const;

  /// Read mesh
  void _read(void);


// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _filename; ///< Name of file

}; // MeshIOSieve

#include "MeshIOSieve.icc" // inline methods

#endif // pylith_meshio_meshiosieve_hh

// End of file 
