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
 * @file pylith/topology/FieldBase.hh
 *
 * @brief Vector field over the vertices or cells of a finite-element
 * mesh or subset of a finite-element mesh.
 *
 * Extends Sieve real general section by adding metadata.
 *
 * We could replace FieldBase, Field, and FieldSubMesh (which
 * implement a field over a finite-element and subset of the
 * finite-element mesh) by templating Field over the mesh, but this
 * would not permit as much insulation against the underlying Sieve
 * implementation. Separation into Field and FieldSubMesh results in
 * complete insulation of objects using Field and FieldSubMesh from
 * the Sieve implementation.
 */

#if !defined(pylith_topology_fieldbase_hh)
#define pylith_topology_fieldbase_hh

// Include directives ---------------------------------------------------
#define NEWPYLITHMESH 1
#include "pylith/utils/sievetypes.hh" // HASA PETSc real_section_type

#include <string> // HASA std::string

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace topology {
    class FieldBase;
    class TestFieldBase;
  } // topology
} // pylith

// FieldBase ----------------------------------------------------------------
class pylith::topology::FieldBase
{ // FieldBase
  friend class TestFieldBase; // unit testing

// PUBLIC ENUMS /////////////////////////////////////////////////////////
public :

  enum VectorFieldEnum {
    SCALAR=0, ///< Scalar.
    VECTOR=1, ///< Vector.
    TENSOR=2, ///< Tensor.
    OTHER=3, ///< Not a scalar, vector, or tensor.
    MULTI_SCALAR=4, ///< Scalar at multiple points.
    MULTI_VECTOR=5, ///< Vector at multiple points.
    MULTI_TENSOR=6, ///< Tensor at multiple points.
    MULTI_OTHER=7, ///< Not a scalar, vector, or tensor at multiple points.
  }; // VectorFieldEnum

  enum DomainEnum {
    VERTICES_FIELD=0, ///< FieldBase over vertices.
    CELLS_FIELD=1, ///< FieldBase over cells.
  }; // omainEnum

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  FieldBase(void);

  /// Destructor.
  ~FieldBase(void);

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

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  double _scale; ///< Dimensional scale associated with field
  std::string _name; ///< Name of field
  VectorFieldEnum _vecFieldType; ///< Type of vector field
  bool _dimensionsOkay; ///< Flag indicating it is okay to dimensionalize

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  FieldBase(const FieldBase&); ///< Not implemented
  const FieldBase& operator=(const FieldBase&); ///< Not implemented

}; // FieldBase

#include "FieldBase.icc"

#endif // pylith_topology_fieldbase_hh


// End of file 
