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
 * @file libsrc/topology/FieldBase.hh
 *
 * @brief Basic information related to a vector field over the
 * vertices or cells of a finite-element mesh.
 */

#if !defined(pylith_topology_fieldbase_hh)
#define pylith_topology_fieldbase_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

// FieldBase ------------------------------------------------------------
class pylith::topology::FieldBase
{ // Field

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
  }; // DomainEnum

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  FieldBase(void); ///< Default constructor.
  ~FieldBase(void); ///< Default destructor.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  FieldBase(const FieldBase&); ///< Not implemented
  const FieldBase& operator=(const FieldBase&); ///< Not implemented

}; // FieldBase

#endif // pylith_topology_fieldbase_hh


// End of file 
