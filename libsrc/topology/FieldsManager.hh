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
 * @file pylith/topology/FieldsManager.hh
 *
 * @brief Object for managing fields associated with the fields
 * defined over a finite-element mesh.
 */

#if !defined(pylith_topology_fieldsmanager_hh)
#define pylith_topology_fieldsmanager_hh

#include <map> // HASA std::map
#include <string> // USES std::string

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

namespace pylith {
  namespace topology {
    class FieldsManager;
    class TestFieldsManager;
  } // topology
} // pylith

class pylith::topology::FieldsManager
{ // FieldsManager
  friend class TestFieldsManager; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  FieldsManager(const ALE::Obj<Mesh>& mesh);

  /// Destructor
  ~FieldsManager(void);

  /** Add field.
   *
   * @param name Name of field.
   */
  void addReal(const char* name);

  /** Get field.
   *
   * @param name Name of field.
   */
  const ALE::Obj<real_section_type>& getReal(const char* name);

  /** Remove field.
   *
   * @param name Name of field.
   */
  void delReal(const char* name);

  /** Copy layout of field to all other fields.
   *
   * @param name Name of field.
   */
  void copyLayout(const char* name);

  /** Copy layout of field to managed fields.
   *
   * @param field Field from which to copy layout.
   */
  void copyLayout(const ALE::Obj<real_section_type>& field);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented
  FieldsManager(const FieldsManager& m);

  /// Not implemented
  const FieldsManager& operator=(const FieldsManager&);

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

  typedef std::map< std::string, ALE::Obj<real_section_type> > map_real_type;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// PETSc mesh associated with fields
  const ALE::Obj<Mesh>& _mesh;

  /// Map for fieldss stored as real fields
  map_real_type _real;

}; // FieldsManager

#endif // pylith_topology_fieldsmanager_hh


// End of file 
