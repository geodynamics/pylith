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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputSolnSubset.hh
 *
 * @brief C++ object for managing output of finite-element data over a
 * subdomain.
 */

#if !defined(pylith_meshio_outputsolnsubset_hh)
#define pylith_meshio_outputsolnsubset_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/topology/Mesh.hh" // ISA OutputManager<Mesh>
#include "pylith/topology/Field.hh" // ISA OutputManager<Field<Mesh>>
#include "OutputManager.hh" // ISA OutputManager

#include <string> // HASA std::string

// OutputSolnSubset -----------------------------------------------------
/** @brief C++ object for managing output of finite-element data over
 * a subdomain.
 */
class pylith::meshio::OutputSolnSubset : public OutputManager
{ // OutputSolnSubset
  friend class TestOutputSolnSubset; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  OutputSolnSubset(void);

  /// Destructor
  ~OutputSolnSubset(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set label identifier for subdomain.
   *
   * @param value Label of subdomain.
   */
  void label(const char* value);

  /** Verify configuration.
   *
   * @param mesh PETSc mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get mesh associated with subdomain.
   *
   * @returns Mesh associated with subdomain.
   */
  const topology::Mesh& subdomainMesh(const topology::Mesh& mesh);
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  OutputSolnSubset(const OutputSolnSubset&); ///< Not implemented.
  const OutputSolnSubset& operator=(const OutputSolnSubset&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _label; ///< Label of subdomain.
  topology::Mesh* _submesh; ///< Mesh of subdomain.

}; // OutputSolnSubset

#endif // pylith_meshio_outputsolnsubset_hh

// End of file 
