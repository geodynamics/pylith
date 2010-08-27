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
 * @brief Set of real fields over a finite-element mesh. All fields are
 * associated with the same points. We use uniform sections so the
 * fiber dimension is set at compile time.
 */

#if !defined(pylith_topology_fields_hh)
#define pylith_topology_fields_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Metadata

#include <string> // USES std::string

// Fields ---------------------------------------------------------------
/// Container for managing multiple fields over a finite-element mesh.
template<typename mesh_type,
	 int fiberDimTotal>
class pylith::topology::FieldsNew
{ // Fields
  friend class TestFieldsNewMesh; // unit testing
  friend class TestFieldsNewSubMesh; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mesh Finite-element mesh.
   */
  Fields(const typename mesh_type& mesh);

  /// Destructor.
  virtual
  ~Fields(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Check if fields contains a given field.
   *
   * @param name Name of field.
   * @return True if fields contains field, false otherwise.
   */
  bool hasField(const char* name) const;

  /** Add field.
   *
   * @param name Name of field.
   * @param label Label for field.
   * @param fiberDim Fiber dimension for field.
   */
  void add(const char* name,
	   const char* label,
	   const int fiberDim,
	   FieldBase::VectorFieldEnum vectorFieldType =FieldBase::OTHER,
	   const double scale =1.0,
	   const bool dimsOkay =false);

  /** Create and allocate Sieve section.
   *
   * @param points Points over which to define section.
   */
  void allocate(const ALE::Obj<label_sequence>& points);

  /** Create and allocate Sieve section.
   *
   * @param points Points over which to define section.
   */
  void allocate(const int_array& points,
		const int fiberDim);

  /** Create and allocate Sieve section.
   *
   * @param domain Type of points over which to define section.
   * @param stratum Stratum depth (for vertices) and height (for cells).
   */
  void allocate(const DomainEnum domain,
		const int stratum =0);

  /** Get field.
   *
   * @param name Name of field.
   */
  const Field<mesh_type>& getField(const char* name) const;
	   
  /** Get field.
   *
   * @param name Name of field.
   */
  Field<mesh_type>& getField(const char* name);
	   
  /** Get mesh associated with fields.
   *
   * @returns Finite-element mesh.
   */
  const typename mesh_type& mesh(void) const;

  /** Return the names of all fields.
   *
   * @returns an array of all field names
   */
  void fieldNames(int *numNames,
		  char ***names);

// PROTECTED STRUCTS ////////////////////////////////////////////////////
protected :

  struct FieldInfo {
    FieldBase::Metadata metadata;
    int fiberDim;
    int fibration; 
  }; // FieldInfo

// PROTECTED TYPEDEFS ///////////////////////////////////////////////////
protected :

  typedef std::map< std::string, FieldInfo > map_type;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  map_type _fields; ///< Fields without constraints over a common set of points.
  ALE::Obj<mesh_type::URealSection> _section; ///< Section containing fields.
  const typename mesh_type& _mesh; ///< Mesh associated with fields.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  FieldsNew(const FieldsNew&); ///< Not implemented
  const FieldsNew& operator=(const FieldsNew&); ///< Not implemented

}; // FieldsNew

#include "FieldsNew.icc"

#endif // pylith_topology_fields_hh


// End of file 
