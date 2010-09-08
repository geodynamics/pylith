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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
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

#define USE_UNIFORMSECTION

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Metadata

#include <string> // USES std::string

// Fields ---------------------------------------------------------------
/// Container for managing multiple fields over a finite-element mesh.
template<typename mesh_type>
class pylith::topology::FieldsNew
{ // Fields
  friend class TestFieldsNewMesh; // unit testing
  friend class TestFieldsNewSubMesh; // unit testing

#if defined(USE_UNIFORMSECTION)
  typedef typename mesh_type::RealUniformSection section_type;
#else
  typedef typename mesh_type::RealSection section_type;
#endif

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mesh Finite-element mesh.
   */
  FieldsNew(const mesh_type& mesh);

  /// Destructor.
  virtual
  ~FieldsNew(void);

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
  void allocate(const ALE::Obj<typename mesh_type::SieveMesh::label_sequence>& points);

  /** Create and allocate Sieve section.
   *
   * @param points Points over which to define section.
   */
  void allocate(const int_array& points);

  /** Create and allocate Sieve section.
   *
   * @param domain Type of points over which to define section.
   * @param stratum Stratum depth (for vertices) and height (for cells).
   */
  void allocate(const FieldBase::DomainEnum domain,
		const int stratum =0);

  /** Get field.
   *
   * @param name Name of field.
   */
  const Field<mesh_type, section_type>& get(const char* name) const;
	   
  /** Get field.
   *
   * @param name Name of field.
   */
  Field<mesh_type, section_type>& get(const char* name);
	   
  /** Get mesh associated with fields.
   *
   * @returns Finite-element mesh.
   */
  const mesh_type& mesh(void) const;

  /** Get section containing fields.
   *
   * @returns Sieve section
   */
  const ALE::Obj<section_type>& section(void) const;

  /** Return the names of all fields.
   *
   * @param numNames Number of fields,
   * @param names Names of fields.
   */
  void fieldNames(int* numNames,
		  std::string** names) const;

// PROTECTED STRUCTS ////////////////////////////////////////////////////
protected :

  struct FieldInfo {
    FieldBase::Metadata metadata;
    int fiberDim;
    int fibration; 
    Field<mesh_type, section_type>* field;
  }; // FieldInfo

// PROTECTED TYPEDEFS ///////////////////////////////////////////////////
protected :

  typedef std::map< std::string, FieldInfo > map_type;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Compute total fiber dimension for section.
   *
   * @returns Fiber dimension.
   */
  int _fiberDim(void) const;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  map_type _fields; ///< Fields without constraints over a common set of points.
  ALE::Obj<section_type> _section; ///< Section containing fields.
  const mesh_type& _mesh; ///< Mesh associated with fields.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  FieldsNew(const FieldsNew&); ///< Not implemented
  const FieldsNew& operator=(const FieldsNew&); ///< Not implemented

}; // FieldsNew

#include "FieldsNew.icc"
#include "FieldsNew.cc"

#endif // pylith_topology_fields_hh


// End of file 
