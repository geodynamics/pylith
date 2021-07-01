// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/PsetFile.hh
 *
 * @brief C++ base class for input/output of LaGriT Pset files.
 */

#if !defined(pylith_meshio_psetfile_hh)
#define pylith_meshio_psetfile_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/array.hh" // HASA int_array
#include <string> // HASA std::string

// PsetFile -------------------------------------------------------------
/// C++ base class for input/output of LaGriT Pset files.
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
