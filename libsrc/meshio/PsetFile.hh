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

#if !defined(pylith_meshio_psetfile_hh)
#define pylith_meshio_psetfile_hh

#include "pylith/utils/array.hh" // USES int_array
#include <string> // HASA std::string

namespace pylith {
  namespace meshio {
  class PsetFile;
  } // meshio
} // pylith

class pylith::meshio::PsetFile
{ // PsetFile

// PUBLIC TYPES /////////////////////////////////////////////////////////
public :

  struct Pset {
    int_array points; ///< Indices of vertices in group
    std::string name; ///< Name of group
    int id; ///< Id of group
  }; // Pset

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Constructor with name of Pset file.
   *
   * @param filename Name of Pset file
   */
  PsetFile(const char* name);

  /// Default destructor.
  ~PsetFile(void);

  /** Is Pset file ascii?
   *
   * @param filename Name of Pset file.
   *
   * @returns True if Pset file is ascii, false otherwise
   */
  static
  bool isAscii(const char* filename);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  std::string _filename; ///< Name of Pset file
  
}; // PsetFile

#endif // pylith_meshio_psetfile_hh


// End of file 
