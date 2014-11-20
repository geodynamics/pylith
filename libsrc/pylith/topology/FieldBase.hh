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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
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

#include <string> // USES std::string

// FieldBase ------------------------------------------------------------
/** @brief Basic information related to a vector field over the
 * vertices or cells of a finite-element mesh.
 */
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
    POINTS_FIELD=2, ///< FieldBase over all points.
    FACES_FIELD=3, ///< FieldBase over faces.
  }; // DomainEnum

// PUBLIC STRUCTS ///////////////////////////////////////////////////////
public :

  struct Metadata {
    std::string label; ///< Label for field.
    VectorFieldEnum vectorFieldType; ///< Type of vector field.
    PylithScalar scale; ///< Dimension scale associated with values.
    bool dimsOkay; ///< Ok to replace nondimensionalized values with dimensionalized values.
  }; // Metadata

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  FieldBase(void); ///< Default constructor.
  ~FieldBase(void); ///< Default destructor.

  /** Get string associated with vector field type.
   *
   * @param value Vector field type.
   * @returns String associated with vector field type.
   */
  static
  const char*
  vectorFieldString(VectorFieldEnum value);

  /** Get string associated with vector field type.
   *
   * @param value String associated with vector field type.
   * @returns Vector field type.
   */
  static
  VectorFieldEnum
  parseVectorFieldString(const char* value);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  FieldBase(const FieldBase&); ///< Not implemented
  const FieldBase& operator=(const FieldBase&); ///< Not implemented

}; // FieldBase

#endif // pylith_topology_fieldbase_hh


// End of file 
