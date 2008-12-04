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
 * @file pylith/topology/Field.hh
 *
 * @brief Vector field over the vertices or cells of a finite-element mesh.
 *
 * Extends Sieve real general section by adding metadata.
 */

#if !defined(pylith_topology_field_hh)
#define pylith_topology_field_hh

// Include directives ---------------------------------------------------
#define NEWPYLITHMESH 1
#include "pylith/utils/sievetypes.hh" // HASA PETSc real_section_type

#include <string> // HASA std::string

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace topology {
    class Field;
    class TestField;
  } // topology
} // pylith

// Field ----------------------------------------------------------------
class pylith::topology::Field
{ // Field
  friend class TestField; // unit testing

// PUBLIC ENUMS /////////////////////////////////////////////////////////
public :

  enum VectorFieldEnum {
    SCALAR, ///< Scalar.
    VECTOR, ///< Vector.
    TENSOR, ///< Tensor.
    OTHER, ///< Not a scalar, vector, or tensor.
    MULTI_SCALAR, ///< Scalar at multiple points.
    MULTI_VECTOR, ///< Vector at multiple points.
    MULTI_TENSOR, ///< Tensor at multiple points.
    MULTI_OTHER, ///< Not a scalar, vector, or tensor at multiple points.
  }; // VectorFieldEnum

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mesh Sieve mesh.
   */
  Field(const ALE::Obj<SieveMesh>& mesh);

  /// Destructor.
  ~Field(void);

  /** Get Sieve section.
   *
   * @returns Sieve section.
   */
  const ALE::Obj<SieveRealSection>& section(void) const;

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

  /** Set spatial dimension of domain.
   *
   * @param value Spatial dimension of domain.
   */
  void spaceDim(const int value);

  /** Get spatial dimension of domain.
   *
   * @returns Spatial dimension of domain.
   */
  int spaceDim(void) const;

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

  /** Create section with same layout (fiber dimension and
   * constraints) as another section. This allows the layout data
   * structures to be reused across multiple fields, reducing memory
   * usage.
   *
   * @param sec Section defining layout.
   */
  void copyLayout(const Field& src);

  /// Clear variables associated with section.
  void clear(void);

  /** Dimensionalize field. Throws runtime_error if field is not
   * allowed to be dimensionalized.
   */
  void dimensionalize(void);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  const ALE::Obj<SieveMesh>& _mesh; ///< Mesh associated with section
  ALE::Obj<SieveRealSection> _section; ///< Real section with data

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  double _scale; ///< Dimensional scale associated with field
  std::string _name; ///< Name of field
  int _spaceDim; ///< Spatial dimension of domain
  VectorFieldEnum _vecFieldType; ///< Type of vector field
  bool _dimensionsOkay; ///< Flag indicating it is okay to dimensionalize

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  Field(const Field&); ///< Not implemented
  const Field& operator=(const Field&); ///< Not implemented

}; // Field

#include "Field.icc"

#endif // pylith_topology_field_hh


// End of file 
