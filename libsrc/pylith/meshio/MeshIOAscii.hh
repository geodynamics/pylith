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

#include "spatialdata/utils/utilsfwd.hh" // USES LineParser

#include <iosfwd> // USES std::istream, std::ostream
#include <string> // HASA std::string

class pylith::meshio::MeshIOAscii : public MeshIO {
    friend class TestMeshIOAscii; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIOAscii(void);

    /// Destructor
    ~MeshIOAscii(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set filename for ASCII file.
     *
     * @param filename Name of file
     */
    void setFilename(const char* name);

    /** Get filename of ASCII file.
     *
     * @returns Name of file
     */
    const char* getFilename(void) const;

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    void _write(void) const;

    /// Read mesh
    void _read(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::string _filename; ///< Name of file
    bool _useIndexZero; ///< Flag indicating if indicates start at 0 (T) or 1 (F)

}; // MeshIOAscii

// End of file
