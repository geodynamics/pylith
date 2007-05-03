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

#if !defined(pylith_meshio_gmvfile_hh)
#define pylith_meshio_gmvfile_hh

#include <string> // HASA std::string

namespace pylith {
  namespace meshio {
  class GMVFile;
  } // meshio
} // pylith

class pylith::meshio::GMVFile
{ // GMVFile

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Constructor with name of GMV file.
   *
   * @param filename Name of GMV file
   */
  GMVFile(const char* name);

  /// Default destructor.
  ~GMVFile(void);

  /** Is GMV file ascii?
   *
   * @param filename Name of GMV file.
   *
   * @returns True if GMV file is ascii, false otherwise
   */
  static
  bool isAscii(const char* filename);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  std::string _filename; ///< Name of GMV file
  
}; // GMVFile

#endif // pylith_meshio_gmvfile_hh


// End of file 
