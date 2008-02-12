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
 * @file pylith/utils/realspaces.hh
 *
 * @brief Categories of real vector space types.
 */

#if !defined(pylith_utils_vectorfields_hh)
#define pylith_utils_vectorfields_hh

namespace pylith {

  /// Enumeration of types of vector fields.
  enum VectorFieldEnum {
    SCALAR_FIELD, ///< Scalar field
    VECTOR_FIELD, ///< Vector field
    TENSOR_FIELD, ///< Tensor field
    OTHER_FIELD ///< Other type of field
  }; // VectorFieldEnum

} // pylith

#endif // pylith_utils_vectorfields_hh


// End of file 
