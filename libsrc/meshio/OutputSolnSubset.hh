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
 * @file pylith/meshio/OutputSolnSubset.hh
 *
 * @brief C++ object for manager output of finite-element data over a
 * subdomain.
 */

#if !defined(pylith_meshio_outputsolnsubset_hh)
#define pylith_meshio_outputsolnsubset_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/topology/SubMesh.hh" // ISA OutputManager<SubMesh>
#include "pylith/topology/Mesh.hh" // ISA OutputManager<Field<Mesh>>
#include "pylith/topology/Field.hh" // ISA OutputManager<Field<Mesh>>
#include "OutputManager.hh" // ISA OutputManager

#include <string> // HASA std::string

// OutputSolnSubset -----------------------------------------------------
class pylith::meshio::OutputSolnSubset : 
  public OutputManager<topology::SubMesh, topology::Field<topology::Mesh> >
{ // OutputSolnSubset
  friend class TestOutputSolnSubset; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  OutputSolnSubset(void);

  /// Destructor
  ~OutputSolnSubset(void);

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
  const topology::SubMesh& subdomainMesh(const topology::Mesh& mesh);
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  OutputSolnSubset(const OutputSolnSubset&); ///< Not implemented.
  const OutputSolnSubset& operator=(const OutputSolnSubset&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _label; ///< Label of subdomain.
  topology::SubMesh* _submesh; ///< Mesh of subdomain.

}; // OutputSolnSubset

#endif // pylith_meshio_outputsolnsubset_hh

// End of file 
