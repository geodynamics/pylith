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

#if !defined(pylith_meshio_meshiolagrit_hh)
#define pylith_meshio_meshiolagrit_hh

#include "MeshIO.hh" // ISA MeshIO

#include <string> // HASA std::string

namespace pylith {
  namespace meshio {
    class MeshIOLagrit;
  } // meshio
} // pylith

class pylith::meshio::MeshIOLagrit : public MeshIO
{ // MeshIOLagrit
// PUBLIC METHODS -------------------------------------------------------
public :

  /// Constructor
  MeshIOLagrit(void);

  /// Destructor
  ~MeshIOLagrit(void);

  /** Set filename for mesh GMV file.
   *
   * @param filename Name of file
   */
  void filenameGmv(const char* name);

  /** Get filename of mesh GMV file.
   *
   * @returns Name of file
   */
  const char* filenameGmv(void) const;

  /** Set filename for PSET mesh file.
   *
   * @param filename Name of file
   */
  void filenamePset(const char* name);

  /** Get filename of PSET mesh file.
   *
   * @returns Name of file
   */
  const char* filenamePset(void) const;

  /** Set flag to write ASCII or binary files.
   *
   * @param flag True if writing ASCII, false if writing binary
   */
  void writerAscii(const bool flag);

  /** Get flag for writing ASCII or binary files.
   *
   * @returns True if writing ASCII, false if writing binary.
   */
  bool writeAscii(void) const;

  /** Set flag to flip endian type when reading/writing from binary files.
   *
   * @param flag True if flipping endian, false otherwise
   */
  void flipEndian(const bool flag);

  /** Get flag for flipping endian type when reading/writing from binary files.
   *
   * @returns True if flipping endian, false othewise.
   */
  bool flipEndian(void) const;

// PROTECTED METHODS ----------------------------------------------------
protected :

  /// Write mesh
  void _write(void) const;

  /// Read mesh
  void _read(void);

// PRIVATE METHODS ------------------------------------------------------
private :

  /** Header in ascii GMV file */
  static const char* GMVASCIIHEADER;

  /** Header in binary GMV file */
  static const char* GMVBINARYHEADER;

  /** Is GMV file an ASCII file?
   *
   * @returns true if GMV file is ASCII, false if binary
   */
  bool _isGmvAscii(void) const;

  /// Read ASCII GMV file.
  void _readGmvAscii(void);

  /// Read binary GMV file.
  void _readGmvBinary(void);

  /// Write ASCII GMV file.
  void _writeGmvAscii(void) const;

  /// Read binary GMV file.
  void _writeGmvBinary(void) const;

  /** Is PSET file an ASCII file?
   *
   * @returns true if PSET file is ASCII, false if binary
   */
  bool _isPsetAscii(void) const;

  /// Read ASCII PSET file.
  void _readPsetAscii(void);

  /// Read binary PSET file.
  void _readPsetBinary(void);

  /// Write ASCII PSET file.
  void _writePsetAscii(void) const;

  /// Write binary PSET file.
  void _writePsetBinary(void) const;

  // PRIVATE MEMBERS ----------------------------------------------------
private :

  std::string _filenameGmv; ///< Name of GMV file
  std::string _filenamePset; ///< Name of PSET file
  bool _writeAscii; ///< True if writing ASCII, false if writing binary
  bool _flipEndian; ///< True if need to change endian when reading/writing

}; // MeshIOLagrit

#include "MeshIOLagrit.icc" // inline methods

#endif // pylith_meshio_meshiolagrit_hh

// End of file 
