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
 * @file libsrc/topology/Fields.hh
 *
 * @brief Object for managing fields over a finite-element mesh.
 */

#if !defined(pylith_topology_fields_hh)
#define pylith_topology_fields_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include <string> // USES std::string

// Fields ---------------------------------------------------------------
template<typename field_type>
class pylith::topology::Fields
{ // Fields
  friend class TestFieldsMesh; // unit testing
  friend class TestFieldsSubMesh; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mesh Finite-element mesh.
   */
  Fields(const typename field_type::Mesh& mesh);

  /// Destructor.
  ~Fields(void);

  /** Add field.
   *
   * @param name Name of field.
   */
  void add(const char* name);

  /** Add field.
   *
   * @param name Name of field.
   * @param domain Type of points over which to define field.
   * @param fiberDim Fiber dimension for field.
   */
  void add(const char* name,
	   const typename field_type::DomainEnum domain,
	   const int fiberDim);

  /** Delete field.
   *
   * @param name Name of field.
   */
  void del(const char* name);

  /** Get field.
   *
   * @param name Name of field.
   */
  const field_type& get(const char* name) const;
	   
  /** Get field.
   *
   * @param name Name of field.
   */
  field_type& get(const char* name);
	   
  /** Copy layout to other fields.
   *
   * @param name Name of field to use as template for layout.
   */
  void copyLayout(const char* name);

  /** Get mesh associated with fields.
   *
   * @returns Finite-element mesh.
   */
  const typename field_type::Mesh& mesh(void) const;

// PROTECTED TYPEDEFS ///////////////////////////////////////////////////
protected :

  typedef std::map< std::string, field_type* > map_type;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  map_type _fields;
  const typename field_type::Mesh& _mesh;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Fields(const Fields&); ///< Not implemented
  const Fields& operator=(const Fields&); ///< Not implemented

}; // Fields

#include "Fields.icc"

#endif // pylith_topology_fields_hh


// End of file 
