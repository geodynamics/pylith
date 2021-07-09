// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/pylith/meshio/Xdmf.hh
 *
 * @brief C++ interface to Python Xdmf object.
 */

#if !defined(pylith_meshio_xdmf_hh)
#define pylith_meshio_xdmf_hh

#include "meshiofwd.hh"

/// C++ interface to Python Xdmf object.
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

#endif // pylith_meshio_xdmf_hh


// End of file
