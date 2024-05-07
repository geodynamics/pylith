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

class pylith::meshio::MeshIOPetsc : public MeshIO {
    friend class TestMeshIOPetsc; // unit testing

    enum Format {
        GMSH=0, // PETSc HDF5 format
        HDF5=1, // PETSc HDF5 format
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIOPetsc(void);

    /// Destructor
    ~MeshIOPetsc(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set filename for ASCII file.
     *
     * @param name Name of file
     */
    void setFilename(const char* name);

    /** Get filename of ASCII file.
     *
     * @returns Name of file
     */
    const char* getFilename(void) const;

    /** Set options prefix for this mesh.
     *
     * @param name Options prefix
     */
    void setPrefix(const char* name);

    /** Get options prefix for this mesh.
     *
     * @returns Options prefix
     */
    const char* getPrefix(void) const;

    /** Set mesh format.
     *
     * @param value Mesh format.
     */
    void setFormat(Format value);

    /** Get mesh format.
     *
     * @returns Mesh format.
     */
    Format getFormat(void) const;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    void _write(void) const;

    /// Read mesh
    void _read(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::string _filename; ///< Name of file.
    std::string _prefix; ///< Options prefix for mesh.
    Format _format; ///< Mesh format.

}; // MeshIOPetsc

// End of file
