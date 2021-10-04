// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/MeshIOPETSc.hh
 *
 * @brief C++ input/output manager for PyLith PETSc gmsh mesh files.
 */

#if !defined(pylith_meshio_meshiopetsc_hh)
#define pylith_meshio_meshiopetsc_hh

#include "MeshIO.hh" // ISA MeshIO

#include <string> // HASA std::string

class pylith::meshio::MeshIOPETSc : public MeshIO {
    friend class TestMeshIOPETSc; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIOPETSc(void);

    /// Destructor
    ~MeshIOPETSc(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set filename for ASCII file.
     *
     * @param name Name of file
     */
    void filename(const char* name);

    /** Get filename of ASCII file.
     *
     * @returns Name of file
     */
    const char* filename(void) const;

    /** Set options prefix for this mesh.
     *
     * @param name Options prefix
     */
    void prefix(const char* name);

    /** Get options prefix for this mesh.
     *
     * @returns Options prefix
     */
    const char* prefix(void) const;


    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    void _write(void) const;

    /// Read mesh
    void _read(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::string _filename; ///< Name of file
    std::string _prefix; ///< Options prefix for mesh

}; // MeshIOPETSc

#include "MeshIOPETSc.icc" // inline methods

#endif // pylith_meshio_meshiopetsc_hh

// End of file
