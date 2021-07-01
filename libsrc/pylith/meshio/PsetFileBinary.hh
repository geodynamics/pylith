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
 * @file libsrc/meshio/PsetFileBinary.hh
 *
 * @brief C++ object for input/output of LaGriT binary Pset files.
 */

#if !defined(pylith_meshio_psetfilebinary_hh)
#define pylith_meshio_psetfilebinary_hh

// Include directives ---------------------------------------------------
#include "PsetFile.hh" // ISA PsetFile

#include <iosfwd>

// PsetFileBinary -------------------------------------------------------
/// C++ object for input/output of LaGriT binary Pset files.
class pylith::meshio::PsetFileBinary : public PsetFile
{ // PsetFileBinary

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Constructor with name of Pset file.
   *
   * @param filename Name of Pset file
   * @param flipEndian Flip endian type when reading/writing.
   * @param ioInt32 True if Pset file uses 64-bit integers.
   * @param isRecordHeader32Bit True if Fortran record header size is 32-bit.
   */
  PsetFileBinary(const char* name,
		 const bool flipEndian,
		 const bool ioInt32,
		 const bool isRecordHeader32Bit);

  /// Default destructor 
  ~PsetFileBinary(void);

  /** Get header.
   *
   * @returns Header that appears in binary Pset file
   */
  static const char* header(void);

  /** Read binary Pset file.
   *
   * @param groups Array of point sets.
   */
  void read(std::vector<Pset>* groups);

  /** Write binary Pset file.
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
  void _readPset32(std::ifstream& fin,
		   Pset* group);

  /** Write point set.
   *
   * @param fout Output file stream.
   * @param group Point set
   */
  void _writePset32(std::ofstream& fout,
		    const Pset& group);

  /** Read point set.
   *
   * @param fin Input file stream.
   * @param group Point set
   */
  void _readPset64(std::ifstream& fin,
		   Pset* group);

  /** Write point set.
   *
   * @param fout Output file stream.
   * @param group Point set
   */
  void _writePset64(std::ofstream& fout,
		    const Pset& group);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :
  
  /** Header in binary Pset file */
  static const char* _HEADER;

  int _recordHeaderSize; ///< Size of Fortran record header in bytes.

  bool _flipEndian; ///< True if need to change endian when reading/writing.
  bool _ioInt32; ///< True if I/O uses pset file uses 32-bit integers.

}; // PsetFileInBinary

#endif // pylith_meshio_psetfilebinary


// End of file 
