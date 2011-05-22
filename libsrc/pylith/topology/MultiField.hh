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
 * @file libsrc/topology/MultiField.hh
 *
 * @brief Vector field over the vertices or cells of a finite-element
 * mesh.
 */

#if !defined(pylith_topology_multifield_hh)
#define pylith_topology_multifield_hh

// Include directives ---------------------------------------------------
#include "FieldBase.hh" // ISA FieldBase

#include "pylith/utils/arrayfwd.hh" // HASA int_array
#include "pylith/utils/petscfwd.h" // HASA PetscVec

#include <petscdmmesh.hh>

#include <map> // USES std::map
#include <string> // USES std::string

// MultiField ----------------------------------------------------------------
/** @brief Vector field over the vertices or cells of a finite-element
 * mesh.
 *
 * Extends Sieve real general section by adding metadata.
 */
template<typename mesh_type,
	 typename section_type>
class pylith::topology::MultiField : public FieldBase
{ // MultiField
  friend class TestMultiFieldMesh; // unit testing
  friend class TestMultiFieldSubMesh; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public:

  // Convenience typedefs
  typedef mesh_type Mesh;

  typedef ALE::ISieveVisitor::RestrictVisitor<section_type> RestrictVisitor;
  typedef ALE::ISieveVisitor::UpdateAddVisitor<section_type> UpdateAddVisitor;
  typedef ALE::ISieveVisitor::UpdateAllVisitor<section_type> UpdateAllVisitor;

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private:

  // Convenience typedefs
  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename SieveMesh::label_sequence label_sequence;
  typedef typename section_type::chart_type chart_type;
  typedef typename mesh_type::SieveMesh::point_type point_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mesh Finite-element mesh.
   */
  MultiField(const mesh_type& mesh);

  /// Destructor.
  ~MultiField(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Get mesh associated with field.
   *
   * @returns Finite-element mesh.
   */
  const mesh_type& mesh(void) const;

  /** Get Sieve section.
   *
   * @returns Sieve section.
   */
  const ALE::Obj<section_type>& section(void) const;

  /** Set label of section.
   *
   * @param value Label of section.
   */
  void label(const char* value);

  /** Get the number of sieve points in the chart.
   *
   * @returns the chart size.
   */
  int chartSize(void) const;

  /** Get the number of degrees of freedom.
   *
   * @returns the number of degrees of freedom.
   */
  int sectionSize(void) const;

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
	   FieldBase::VectorFieldEnum vectorFieldType =FieldBase::OTHER,
	   const double scale =1.0,
	   const bool dimsOkay =false);

  /** Get field.
   *
   * @param name Name of field.
   * @returns Field.
   */
  Field<mesh_type>& get(const char* name);
	   
  /** Get index of field in collection of fields.
   *
   * @param name Name of field.
   * @returns Index of field in collection of fields.
   */
  int fieldIndex(const char* name) const;

  /** Get index of first value of field in field.
   *
   * @param fieldIndex Index of field in collection.
   * @param point Point in finite-element mesh.
   * @returns Index of first value of field in section.
   */
  int fieldStartIndex(const int fieldIndex,
		      const point_type point) const;

  /** Get fiber dimension of field in section.
   *
   * @param fieldIndex Index of field in collection.
   * @param point Point in finite-element mesh.
   * @returns Fiber dimension of field in section.
   */
  int fieldFiberDim(const int fieldIndex,
		    const point_type point) const;

  /** Compute total fiber dimension for section.
   *
   * @param point Point in finite-element mesh.
   * @returns Fiber dimension.
   */
  int fiberDim(const point_type point) const;

  /// Clear variables associated with section.
  void clear(void);

  /// Allocate field.
  void allocate(void);

  /** Create section with same layout (fiber dimension and
   * constraints) as another section. This allows the layout data
   * structures to be reused across multiple fields, reducing memory
   * usage.
   *
   * @param sec MultiField defining layout.
   */
  void cloneSection(const MultiField& src);

  /// Complete section by assembling across processors.
  void complete(void);

  /// Zero section values (does not zero constrained values).
  void zero(void);

  /// Zero section values (including constrained values).
  void zeroAll(void);

  /** Copy field values and metadata.
   *
   * @param field MultiField to copy.
   */
  void copy(const MultiField& field);

  /** Copy field values.
   *
   * @param field MultiField to copy.
   */
  void copy(const ALE::Obj<section_type>& field);

  /** Add two fields, storing the result in one of the fields.
   *
   * @param field MultiField to add.
   */
  MultiField& operator+=(const MultiField& field);

  /** Dimensionalize field. Throws runtime_error if field is not
   * allowed to be dimensionalized.
   */
  void dimensionalize(void) const;

  /** Print field to standard out.
   *
   * @param label Label for output.
   */
  void view(const char* label) const;

  /** Create PETSc vector scatter for field. This is used to transfer
   * information from the "global" PETSc vector view to the "local"
   * Sieve section view. The PETSc vector does not contain constrained
   * DOF. Use createScatterWithBC() to include the constrained DOF in
   * the PETSc vector.
   *
   * @param context Label for context associated with vector.
   */
  void createScatter(const char* context ="");

  /** Create PETSc vector scatter for field. This is used to transfer
   * information from the "global" PETSc vector view to the "local"
   * Sieve section view. The PETSc vector does not contain constrained
   * DOF. Use createScatterWithBC() to include the constrained DOF in
   * the PETSc vector.
   *
   * @param numbering Numbering used to select points in section.
   * @param context Label for context associated with vector.
   */
  void createScatter(const typename ALE::Obj<typename SieveMesh::numbering_type> numbering,
		     const char* context ="");

  /** Create PETSc vector scatter for field. This is used to transfer
   * information from the "global" PETSc vector view to the "local"
   * Sieve section view. The PETSc vector includes constrained
   * DOF. Use createScatter() if constrained DOF should be omitted
   * from the PETSc vector.
   *
   * @param context Label for context associated with vector.
   */
  void createScatterWithBC(const char* context ="");

  /** Create PETSc vector scatter for field. This is used to transfer
   * information from the "global" PETSc vector view to the "local"
   * Sieve section view. The PETSc vector includes constrained
   * DOF. Use createScatter() if constrained DOF should be omitted
   * from the PETSc vector.
   *
   * @param numbering Numbering used to select points in section.
   * @param context Label for context associated with vector.
   */
  void createScatterWithBC(const typename ALE::Obj<typename SieveMesh::numbering_type> numbering,
		     const char* context ="");

  /** Get PETSc vector associated with field.
   *
   * @param context Label for context associated with vector.
   * @returns PETSc vector.
   */
  PetscVec vector(const char* context ="");

  /** Get PETSc vector associated with field.
   *
   * @param context Label for context associated with vector.
   * @returns PETSc vector.
   */
  const PetscVec vector(const char* context ="") const;

  /// Scatter section information across processors to update the
  /// PETSc vector view of the field.
  void scatterSectionToVector(const char* context ="") const;

  /** Scatter section information across processors to update the
   * PETSc vector view of the field.
   *
   * @param vector PETSc vector to update.
   */
  void scatterSectionToVector(const PetscVec vector,
			      const char* context ="") const;

  /// Scatter PETSc vector information across processors to update the
  /// Sieve section view of the field.
  void scatterVectorToSection(const char* context ="") const;

  /** Scatter section information across processors to update the
   * PETSc vector view of the field.
   *
   * @param vector PETSc vector used in update.
   */
  void scatterVectorToSection(const PetscVec vector,
			      const char* context ="") const;

  /// Setup split field with all entries set to a default space of 0.
  void splitDefault(void);

// PRIVATE STRUCTS //////////////////////////////////////////////////////
private :

  struct FieldInfo {
    FieldBase::Metadata metadata; ///< Metadata for field.
    int fieldIndex;  ///< Index associated with field.
    Field<mesh_type>* field; ///< Single field.
  }; // FieldInfo

  /// Data structures used in scattering to/from PETSc Vecs.
  struct ScatterInfo {
    PetscVec vector; ///< PETSc vector associated with field.
    PetscVecScatter scatter; ///< PETSc scatter associated with field.
    PetscVec scatterVec; ///< PETSC vector used in scattering.
  }; // ScatterInfo

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

  typedef std::map<std::string, FieldInfo> metadata_map_type;
  typedef std::map<std::string, ScatterInfo> scatter_map_type;


// PRIVATE METHODS //////////////////////////////////////////////////////
private :

  /** Get scatter for given context.
   *
   * @param context Context for scatter.
   * @param createOk If true, okay to create new scatter for
   * context, if false will throw std::runtime_error if scatter for
   * context doesn't already exist.
   */
  ScatterInfo& _getScatter(const char* context,
			   const bool createOk =false);

  /** Get scatter for given context.
   *
   * @param context Context for scatter.
   */
  const ScatterInfo& _getScatter(const char* context) const;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  const mesh_type& _mesh; ///< Mesh associated with section.
  ALE::Obj<section_type> _section; ///< Real section with data.
  metadata_map_type _fields; ///< Metadata for fields in section.
  scatter_map_type _scatters; ///< Collection of scatters.


// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  MultiField(const MultiField&); ///< Not implemented
  const MultiField& operator=(const MultiField&); ///< Not implemented

}; // MultiField

#include "MultiField.icc"
#include "MultiField.cc"

#endif // pylith_topology_multifield_hh


// End of file 
