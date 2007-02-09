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
 *
 * :QUESTION: Should there be only one parameter manager for all
 * element families or one per element family? If we can easily create
 * a section that pertains only to an element family, then one per
 * family would make the most sense.
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
{ // Integrator
  friend class TestParameterManager; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

  typedef ALE::Mesh Mesh;
  typedef Mesh::topology_type topology_type;
  typedef topology_type::point_type point_type;
  typedef Mesh::real_section_type real_section_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  ParameterManager(void);

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

  /// Map for parameters stored as real fields
  map_real_type _real;

}; // ParameterManager

#endif // pylith_feassemble_parametermanager_hh


// End of file 
