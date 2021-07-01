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
 * @file libsrc/meshio/MeshIOLagrit.hh
 *
 * @brief Input/output manager for LaGriT GMV and Pset files.
 */

#if !defined(pylith_meshio_meshiolagrit_hh)
#define pylith_meshio_meshiolagrit_hh

#include "MeshIO.hh" // ISA MeshIO

#include <string> // HASA std::string

class pylith::meshio::MeshIOLagrit : public MeshIO {
    friend class TestMeshIOLagrit; // unit testing

    // PUBLIC METHODS
    // //////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIOLagrit(void);

    /// Destructor
    ~MeshIOLagrit(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

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

    /** Set flag indicating LaGriT Pset files use 32-bit integers.
     *
     * @param flag True if using 32-bit integers, false if using 64-bit integers.
     */
    void ioInt32(const bool flag);

    /** Get flag indicating LaGriT Pset files use 32-bit integers.
     *
     * @returns True if using 32-bit integers, false if using 64-bit integers.
     */
    bool ioInt32(void) const;

    /** Set Fortran record header size flag.
     *
     * @param flag True if Fortran record header size is 32-bit, false if 64-bit.
     */
    void isRecordHeader32Bit(const bool flag);

    /** Get Fortran record header size flag.
     *
     * @param returns True if Fortran record header size is 32-bit,
     *   false if 64-bit.
     */
    bool isRecordHeader32Bit(void) const;

    // PROTECTED METHODS
    // ///////////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    void _write(void) const;

    /// Read mesh
    void _read(void);

    // PRIVATE METHODS
    // /////////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /** Reorder vertices in cells from ASCII GMV files to match PyLith
     * conventions.
     *
     * @param cells Array of vertex indices for each cell [numCells*numCorners].
     * @param numCells Number of cells.
     * @param numCorners Number of vertices per cell.
     * @param meshDim Spatial dimension of mesh.
     */
    static
    void _orientCellsAscii(int_array* const cells,
                           const int numCells,
                           const int numCorners,
                           const int meshDim);

    /** Reorder vertices in cells from binary GMV files to match PyLith
     * conventions.
     *
     * @param cells Array of vertex indices for each cell [numCells*numCorners].
     * @param numCells Number of cells.
     * @param numCorners Number of vertices per cell.
     * @param meshDim Spatial dimension of mesh.
     */
    static
    void _orientCellsBinary(int_array* const cells,
                            const int numCells,
                            const int numCorners,
                            const int meshDim);

    // PRIVATE MEMBERS
    // /////////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::string _filenameGmv; ///< Name of GMV file.
    std::string _filenamePset; ///< Name of PSET file.
    bool _writeAscii; ///< True if writing ASCII, false if writing binary.
    bool _flipEndian; ///< True if need to change endian when reading/writing.
    bool _ioInt32; ///< True if using 64-bit integers in Pset files.
    bool _isRecordHeader32Bit; ///< True if Fortran record header is 32-bit.

}; // MeshIOLagrit

#include "MeshIOLagrit.icc" // inline methods

#endif // pylith_meshio_meshiolagrit_hh

// End of file
