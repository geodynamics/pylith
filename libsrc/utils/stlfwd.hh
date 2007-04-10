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
 * @file pylith/utils/stlfwd.hh
 *
 * @brief Forward declarations for C++ STL objects.
 *
 * For simple types (i.e., int and double) std::valarray provides some
 * features that std::vector does not have, such as operating on the
 * whole array at once.
 */

#if !defined(pylith_utils_stlfwd_hh)
#define pylith_utils_stlfwd_hh

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

  /// Alias for std::valarray<int>
  typedef std::valarray<int> int_array;

  /// Alias for std::valarray<double>
  typedef std::valarray<double> double_array;

} // pylith


#endif // pylith_utils_stlfwd_hh

// End of file
