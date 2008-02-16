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

#include "OutputManager.hh" // ISA OutputManager

#include "pylith/utils/sievetypes.hh" // HASA ALE::Mesh
#include <string> // HASA std::string

namespace pylith {
  namespace meshio {
    class OutputSolnSubset;
  } // meshio
} // pylith

class pylith::meshio::OutputSolnSubset : public OutputManager
{ // OutputSolnSubset

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

  /** Get mesh associated with subdomain.
   *
   * @returns Mesh associated with subdomain.
   */
  const ALE::Obj<Mesh>& subdomainMesh(const ALE::Obj<Mesh>& mesh);
  

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  OutputSolnSubset(const OutputSolnSubset&); ///< Not implemented.
  const OutputSolnSubset& operator=(const OutputSolnSubset&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  std::string _label; ///< Label of subdomain.
  ALE::Obj<Mesh> _mesh; ///< Mesh of subdomain.

}; // OutputSolnSubset

#endif // pylith_meshio_outputsolnsubset_hh

// End of file 
