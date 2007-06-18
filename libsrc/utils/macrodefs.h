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
 * @file pylith/utils/macrodefs.hh
 *
 * @brief Macro definitions for PyLith.
 */

#if !defined(pylith_utils_macrodefs_h)
#define pylith_utils_macrodefs_h

#if !defined(CALL_MEMBER_FN)
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
#endif

#endif // pylith_utils_macro_defs_h


// End of file
