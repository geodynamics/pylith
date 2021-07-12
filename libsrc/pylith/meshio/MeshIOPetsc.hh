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
 * @file libsrc/meshio/MeshIOPetsc.hh
 *
 * @brief C++ input/output manager for PyLith PETSc gmsh mesh files.
 */

#if !defined(pylith_meshio_meshiopetsc_hh)
#define pylith_meshio_meshiopetsc_hh

#include "MeshIO.hh" // ISA MeshIO

#include <string> // HASA std::string

class pylith::meshio::MeshIOPetsc : public MeshIO {
    friend class TestMeshIOPetsc; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIOPetsc(void);

    /// Destructor
    ~MeshIOPetsc(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);


    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    void _write(void) const;

    /// Read mesh
    void _read(void);

    /// Read mesh
    void _label(void);
    
    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    bool _useIndexZero; ///< Flag indicating if indicates start at 0 (T) or 1 (F)

}; // MeshIOPetsc


#endif // pylith_meshio_meshiopetsc_hh

// End of file
