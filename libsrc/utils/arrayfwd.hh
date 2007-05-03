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
 * @file pylith/utils/arrayfwd.hh
 *
 * @brief Forward declarations for PyLith array objects.
 *
 * These are generally just forward declarations for C++ STL objects.
 *
 * For simple types (i.e., int and double) std::valarray provides some
 * features that std::vector does not have, such as operating on the
 * whole array at once.
 */

#if !defined(pylith_utils_arrayfwd_hh)
#define pylith_utils_arrayfwd_hh

#include <string> // USES std::string

/// Forward declaration of STL vector
namespace std {
  // std::vector
  template<typename T> class allocator;
  template<typename T, typename U> class vector;

  // std::valarray
  template<typename T> class valarray;
} // std

/// Aliases 
namespace pylith {
  /// Alias for std::vector<int>
  typedef std::vector<int, std::allocator<int> > int_vector;

  /// Alias for std::vector<double>
  typedef std::vector<double, std::allocator<double> > double_vector;

  /// Alias for std::vector<std::string>
  typedef std::vector<std::string, std::allocator<std::string> > string_vector;

  /// Alias for std::valarray<int>
  typedef std::valarray<int> int_array;

  /// Alias for std::valarray<float>
  typedef std::valarray<float> float_array;

  /// Alias for std::valarray<double>
  typedef std::valarray<double> double_array;

} // pylith


#endif // pylith_utils_arrayfwd_hh

// End of file
