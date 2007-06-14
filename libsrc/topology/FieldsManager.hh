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
#include <vector> // USES std::vector

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

  /** Set fiber dimension for field.
   *
   * @param name Name of field.
   * @param fiberDim Fiber dimension
   * @param points Name of set of points in mesh over which section is defined
   *   { "vertices", "cells" }
   */
  void setFiberDimension(const char* name,
			 const int fiberDim,
			 const char* points);

  /** Allocate field.
   *
   * @param name Name of field.
   */
  void allocate(const char* name);

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

  /** Set name of solution field.
   *
   * @param name Name of field that is the solution.
   */
  void solutionField(const char* name);

  /** Get solution field.
   *
   * @returns Solution field.
   */
  const ALE::Obj<real_section_type>& getSolution(void);

  /** Create history manager for a subset of the managed fields.
   *
   * @param fields Fields in history (first is most recent).
   * @param size Number of fields in history.
   */
  void createHistory(const char** fields,
		     const int size);

  /** Shift fields in history. Handles to fields are shifted so that
   *  the most recent values become associated with the second most
   *  recent item in the history, etc.
   */
  void shiftHistory(void);

  /** Get field in history by position.
   *
   * @param index Index in history [0=most recent, 1=previous, etc]
   * @returns Field in history.
   */
  const ALE::Obj<real_section_type>& getHistoryItem(const int index);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented
  FieldsManager(const FieldsManager&);

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

  /// Name of field that corresponds to the solution.
  std::string _solutionName;

  /// History manager for a subset of the fields
  std::vector<std::string> _history;

}; // FieldsManager

#endif // pylith_topology_fieldsmanager_hh


// End of file 
