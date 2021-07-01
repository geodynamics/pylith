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
 * @file libsrc/meshio/PsetFileAscii.hh
 *
 * @brief C++ object for input/output of LaGriT binary Pset files.
 */

#if !defined(pylith_meshio_psetfileascii_hh)
#define pylith_meshio_psetfileascii_hh

// Include directives ---------------------------------------------------
#include "PsetFile.hh" // ISA PsetFile

#include <iosfwd>

// PsetFileAscii --------------------------------------------------------
/// C++ object for input/output of LaGriT ascii Pset files.
class pylith::meshio::PsetFileAscii : public PsetFile
{ // PsetFileAscii

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Constructor with name of Pset file.
   *
   * @param filename Name of Pset file
   */
  PsetFileAscii(const char* name);

  /// Default destructor 
  ~PsetFileAscii(void);

  /** Get header.
   *
   * @returns Header that appears in ASCII Pset file
   */
  static const char* header(void);

  /** Read ASCII Pset file.
   *
   * @param groups Array of point sets.
   */
  void read(std::vector<Pset>* groups);

  /** Write ASCII Pset file.
   *
   * @param groups Array of point sets.
   */
  void write(const std::vector<Pset>& groups);

// PRIVATE METHODS //////////////////////////////////////////////////////
private :
  
  /** Read header.
   *
   * @param fin Input file stream
   */
  void _readHeader(std::ifstream& fin);

  /** Write header.
   *
   * @param fout Output file stream
   */
  void _writeHeader(std::ofstream& fout);

  /** Read point set.
   *
   * @param fin Input file stream.
   * @param group Point set
   */
  void _readPset(std::ifstream& fin,
		 Pset* group);

  /** Write point set.
   *
   * @param fout Output file stream.
   * @param group Point set
   */
  void _writePset(std::ofstream& fout,
		  const Pset& group);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :
  
  /** Header in ascii Pset file */
  static const char* _HEADER;

}; // PsetFileInAscii

#include "PsetFileAscii.icc" // inline methods

#endif // pylith_meshio_psetfileascii


// End of file 
