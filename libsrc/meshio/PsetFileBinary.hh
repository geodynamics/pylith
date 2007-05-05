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

#if !defined(pylith_meshio_psetfilebinary_hh)
#define pylith_meshio_psetfilebinary_hh

#include "PsetFile.hh" // ISA PsetFile

#include "pylith/utils/array.hh" // USES int_array, string_array, std::vector
#include <iosfwd>

namespace pylith {
  namespace meshio {
    class PsetFileBinary;
  } // meshio
} // pylith

class pylith::meshio::PsetFileBinary : public PsetFile
{ // PsetFileBinary

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /** Constructor with name of Pset file.
   *
   * @param filename Name of Pset file
   * @param flipEndian Flip endian type when reading/writing.
   */
  PsetFileBinary(const char* name,
		 const bool flipEndian);

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
  
  /** Header in binary Pset file */
  static const char* _HEADER;

  bool _flipEndian; ///< True if need to change endian when reading/writing

}; // PsetFileInBinary

#endif // pylith_meshio_psetfilebinary


// End of file 
