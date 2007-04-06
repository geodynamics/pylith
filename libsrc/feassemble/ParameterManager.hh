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
 * @file pylith/feassemble/ParameterManager.hh
 *
 * @brief Object for managing fields associated with parameters for
 * finite-elements.
 *
 * The parameter manager stores the fields associated with parameters
 * for physical properties or boundary conditions.
 */

#if !defined(pylith_feassemble_parametermanager_hh)
#define pylith_feassemble_parametermanager_hh

#include <map> // HASA std::map
#include <string> // USES std::string

#include <petscmesh.h> // USES Mesh

namespace pylith {
  namespace feassemble {
    class ParameterManager;
    class TestParameterManager;
  } // feassemble
} // pylith

class pylith::feassemble::ParameterManager
{ // ParameterManager
  friend class TestParameterManager; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef ALE::Mesh               Mesh;
  typedef Mesh::point_type        point_type;
  typedef Mesh::real_section_type real_section_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ParameterManager(const ALE::Obj<Mesh>& mesh);

  /// Destructor
  ~ParameterManager(void);

  /** Add parameter.
   *
   * @param name Name of parameter field
   */
  void addReal(const char* name);

  /** Get parameter.
   *
   * @param name Name of parameter field
   */
  const ALE::Obj<real_section_type>& getReal(const char* name);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /// Not implemented
  ParameterManager(const ParameterManager& m);

  /// Not implemented
  const ParameterManager& operator=(const ParameterManager&);

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

  typedef std::map< std::string, ALE::Obj<real_section_type> > map_real_type;


// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// PETSc mesh associated with fields
  const ALE::Obj<Mesh>& _mesh;

  /// Map for parameters stored as real fields
  map_real_type _real;

}; // ParameterManager

#endif // pylith_feassemble_parametermanager_hh


// End of file 
