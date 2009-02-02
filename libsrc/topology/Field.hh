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

#include "FieldBase.hh" // ISA FieldBase

#if !defined(pylith_topology_field_hh)
#define pylith_topology_field_hh

// Include directives ---------------------------------------------------
#define NEWPYLITHMESH 1
#include "pylith/utils/sievetypes.hh" // HASA PETSc real_section_type

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace topology {
    class Field;
    class TestField;

    class Mesh; // HASA Mesh
  } // topology
} // pylith

// Field ----------------------------------------------------------------
class pylith::topology::Field : public FieldBase
{ // Field
  friend class TestField; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mesh Finite-element mesh.
   */
  Field(const Mesh& mesh);

  /// Destructor.
  ~Field(void);

  /** Get Sieve section.
   *
   * @returns Sieve section.
   */
  const ALE::Obj<MeshRealSection>& section(void) const;

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
  void newSection(const ALE::Obj<SieveMesh::label_sequence>& points,
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
  void newSection(const MeshRealSection::chart_type& chart,
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

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  const Mesh& _mesh; ///< Mesh associated with section
  ALE::Obj<MeshRealSection> _section; ///< Real section with data

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Field(const Field&); ///< Not implemented
  const Field& operator=(const Field&); ///< Not implemented

}; // Field

#include "Field.icc"

#endif // pylith_topology_field_hh


// End of file 
