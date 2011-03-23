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
 * @file libsrc/topology/Field.hh
 *
 * @brief Vector field over the vertices or cells of a finite-element
 * mesh.
 */

#if !defined(pylith_topology_field_hh)
#define pylith_topology_field_hh

// Include directives ---------------------------------------------------
#include "FieldBase.hh" // ISA FieldBase

#include "pylith/utils/arrayfwd.hh" // HASA int_array
#include "pylith/utils/petscfwd.h" // HASA PetscVec

#include <petscdmmesh.hh>

// Field ----------------------------------------------------------------
/** @brief Vector field over the vertices or cells of a finite-element
 * mesh.
 *
 * Extends Sieve real general section by adding metadata.
 */
template<typename mesh_type>
class pylith::topology::Field : public FieldBase
{ // Field
  friend class TestFieldMesh; // unit testing
  friend class TestFieldSubMesh; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public:

  // Convenience typedefs
  typedef mesh_type Mesh;

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private:

  // Convenience typedefs
  typedef typename mesh_type::RealSection RealSection;
  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename SieveMesh::label_sequence label_sequence;
  typedef typename RealSection::chart_type chart_type;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mesh Finite-element mesh.
   */
  Field(const mesh_type& mesh);

  /// Destructor.
  ~Field(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Get Sieve section.
   *
   * @returns Sieve section.
   */
  const ALE::Obj<RealSection>& section(void) const;

  /** Get mesh associated with field.
   *
   * @returns Finite-element mesh.
   */
  const mesh_type& mesh(void) const;

  /** Set label for field.
   *
   * @param value Label for field.
   */
  void label(const char* value);

  /** Get label for field.
   *
   * @returns Label for field.
   */
  const char* label(void) const;

  /** Set vector field type
   *
   * @param value Type of vector field.
   */
  void vectorFieldType(const VectorFieldEnum value);

  /** Get vector field type
   *
   * @returns Type of vector field.
   */
  VectorFieldEnum vectorFieldType(void) const;

  /** Set scale for dimensionalizing field.
   *
   * @param value Scale associated with field.
   */
  void scale(const double value);

  /** Get scale for dimensionalizing field.
   *
   * @returns Scale associated with field.
   */
  double scale(void) const;

  /** Set flag indicating whether it is okay to dimensionalize field.
   *
   * @param value True if it is okay to dimensionalize field.
   */
  void addDimensionOkay(const bool value);

  /** Set flag indicating whether it is okay to dimensionalize field.
   *
   * @param value True if it is okay to dimensionalize field.
   */
  bool addDimensionOkay(void) const;

  /** Get spatial dimension of domain.
   *
   * @returns Spatial dimension of domain.
   */
  int spaceDim(void) const;

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

  /// Create sieve section.
  void newSection(void);

  /** Create sieve section and set chart and fiber dimesion for
   * sequence of points.
   *
   * @param points Points over which to define section.
   * @param dim Fiber dimension for section.
   */
  void newSection(const ALE::Obj<label_sequence>& points,
		  const int fiberDim);

  /** Create sieve section and set chart and fiber dimesion for a list
   * of points.
   *
   * @param points Points over which to define section.
   * @param dim Fiber dimension for section.
   */
  void newSection(const int_array& points,
		  const int fiberDim);

  /** Create sieve section and set chart and fiber dimesion.
   *
   * @param domain Type of points over which to define section.
   * @param dim Fiber dimension for section.
   * @param stratum Stratum depth (for vertices) and height (for cells).
   */
  void newSection(const DomainEnum domain,
		  const int fiberDim,
		  const int stratum =0);

  /** Create section using src field as template with given fiber dimension.
   *
   * @param sec Field defining layout.
   * @param fiberDim Fiber dimension.
   */
  void newSection(const Field& src,
		  const int fiberDim);

  /** Create section with same layout (fiber dimension and
   * constraints) as another section. This allows the layout data
   * structures to be reused across multiple fields, reducing memory
   * usage.
   *
   * @param sec Field defining layout.
   */
  void cloneSection(const Field& src);

  /// Clear variables associated with section.
  void clear(void);

  /// Allocate field.
  void allocate(void);

  /// Zero section values (does not zero constrained values).
  void zero(void);

  /// Zero section values (including constrained values).
  void zeroAll(void);

  /// Complete section by assembling across processors.
  void complete(void);

  /** Copy field values and metadata.
   *
   * @param field Field to copy.
   */
  void copy(const Field& field);

  /** Copy field values.
   *
   * @param field Field to copy.
   */
  void copy(const ALE::Obj<typename mesh_type::RealSection>& field);

  /** Add two fields, storing the result in one of the fields.
   *
   * @param field Field to add.
   */
  Field& operator+=(const Field& field);

  /** Dimensionalize field. Throws runtime_error if field is not
   * allowed to be dimensionalized.
   */
  void dimensionalize(void) const;

  /** Print field to standard out.
   *
   * @param label Label for output.
   */
  void view(const char* label) const;

  /// Create PETSc vector for field.
  void createVector(void);

  /** Get PETSc vector associated with field.
   *
   * @returns PETSc vector.
   */
  PetscVec vector(void);

  /** Get PETSc vector associated with field.
   *
   * @returns PETSc vector.
   */
  const PetscVec vector(void) const;

  /// Create PETSc vector scatter for field. This is used to transfer
  /// information from the "global" PETSc vector view to the "local"
  /// Sieve section view.
  void createScatter(void);

  /// Scatter section information across processors to update the
  /// PETSc vector view of the field.
  void scatterSectionToVector(void) const;

  /** Scatter section information across processors to update the
   * PETSc vector view of the field.
   *
   * @param vector PETSc vector to update.
   */
  void scatterSectionToVector(const PetscVec vector) const;

  /// Scatter PETSc vector information across processors to update the
  /// Sieve section view of the field.
  void scatterVectorToSection(void) const;

  /** Scatter section information across processors to update the
   * PETSc vector view of the field.
   *
   * @param vector PETSc vector used in update.
   */
  void scatterVectorToSection(const PetscVec vector) const;

  /// Setup split field with all entries set to a default space of 0.
  void splitDefault(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _scale; ///< Dimensional scale associated with field.
  std::string _label; ///< Label for field.
  const mesh_type& _mesh; ///< Mesh associated with section.
  ALE::Obj<RealSection> _section; ///< Real section with data.
  PetscVec _vector; ///< PETSc vector associated with field.
  PetscVecScatter _scatter; ///< PETSc scatter associated with field.
  PetscVec _scatterVec; ///< PETSC vector used in scattering.
  VectorFieldEnum _vecFieldType; ///< Type of vector field.
  bool _dimensionsOkay; ///< Flag indicating it is okay to dimensionalize.


// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Field(const Field&); ///< Not implemented
  const Field& operator=(const Field&); ///< Not implemented

}; // Field

#include "Field.icc"
#include "Field.cc"

#endif // pylith_topology_field_hh


// End of file 
