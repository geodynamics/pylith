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
 * @file libsrc/topology/Field.hh
 *
 * @brief Vector field over the vertices or cells of a finite-element
 * mesh.
 *
 * Extends Sieve real general section by adding metadata.
 */

#if !defined(pylith_topology_field_hh)
#define pylith_topology_field_hh

// Include directives ---------------------------------------------------
#include "FieldBase.hh" // ISA FieldBase

#include "pylith/utils/petscfwd.h" // HASA PetscVec

#include <petscmesh.hh>

// Field ----------------------------------------------------------------
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

  /** Set name of field.
   *
   * @param value Name of field.
   */
  void name(const char* value);

  /** Get name of field.
   *
   * @returns Name of field.
   */
  const char* name(void) const;

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

  /// Create sieve section.
  void newSection(void);

  /** Create sieve section and set chart and fiber dimesion.
   *
   * @param points Points over which to define section.
   * @param dim Fiber dimension for section.
   */
  void newSection(const ALE::Obj<label_sequence>& points,
		  const int fiberDim);

  /** Create sieve section and set chart and fiber dimesion.
   *
   * @param domain Type of points over which to define section.
   * @param dim Fiber dimension for section.
   */
  void newSection(const DomainEnum domain,
		  const int fiberDim);

  /** Create section given chart. This allows a chart to be reused
   * across multiple fields, reducing memory usage.
   *
   * @param chart Chart defining points over which section is defined.
   * @param fiberDim Fiber dimension.
   */
  void newSection(const chart_type& chart,
		  const int fiberDim);

  /** Create section with same layout (fiber dimension and
   * constraints) as another section. This allows the layout data
   * structures to be reused across multiple fields, reducing memory
   * usage.
   *
   * @param sec Section defining layout.
   */
  void newSection(const Field& src);

  /// Clear variables associated with section.
  void clear(void);

  /// Allocate field.
  void allocate(void);

  /// Zero section values.
  void zero(void);

  /// Complete section by assembling across processors.
  void complete(void);

  /** Copy field values and metadata.
   *
   * @param field Field to copy.
   */
  void copy(const Field& field);

  /** Add two fields, storing the result in one of the fields.
   *
   * @param field Field to add.
   */
  void operator+=(const Field& field);

  /** Dimensionalize field. Throws runtime_error if field is not
   * allowed to be dimensionalized.
   */
  void dimensionalize(void);

  /** Print field to standard out.
   *
   * @param label Label for output.
   */
  void view(const char* label);

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

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _scale; ///< Dimensional scale associated with field
  std::string _name; ///< Name of field
  const mesh_type& _mesh; ///< Mesh associated with section
  ALE::Obj<RealSection> _section; ///< Real section with data
  PetscVec _vector;
  VectorFieldEnum _vecFieldType; ///< Type of vector field
  bool _dimensionsOkay; ///< Flag indicating it is okay to dimensionalize


// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Field(const Field&); ///< Not implemented
  const Field& operator=(const Field&); ///< Not implemented

}; // Field

#include "Field.icc"
#include "Field.cc"

#endif // pylith_topology_field_hh


// End of file 
