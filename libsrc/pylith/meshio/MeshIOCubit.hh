// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/meshio/MeshIO.hh" // ISA MeshIO

#include <string> // HASA std::string

class pylith::meshio::MeshIOCubit : public MeshIO {
    friend class TestMeshIOCubit; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIOCubit(void);

    /// Destructor
    ~MeshIOCubit(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set filename for Cubit file.
     *
     * @param filename Name of file
     */
    void setFilename(const char* name);

    /** Get filename of Cubit file.
     *
     * @returns Name of file
     */
    const char* getFilename(void) const;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    void _write(void) const;

    /// Read mesh
    void _read(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Read groups of vertices.
     *
     * @param[inout] fileIn Input file.
     */
    void _readNodeSets(pylith::meshio::ExodusII& fileIn);

    /** Read groups of faces.
     *
     * @param[inout] fileIn Input file.
     */
    void _readSideSets(pylith::meshio::ExodusII& fileIn);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::string _filename; ///< Name of file

}; // MeshIOCubit

// End of file
