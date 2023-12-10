// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/meshio/meshiofwd.hh"

class pylith::meshio::Xdmf {
    friend class TestXdmf; // Unit testing

public:

    // PUBLIC METHODS /////////////////////////////////////////////////////

    /** Use Python Xdmf object to write Xdmf file corresponding to HDF5 file.
     *
     * @param[in] filenameH5 Name of HDF5 file.
     */
    static
    void write(const char* filenameH5);

}; // class Xdmf

// End of file
