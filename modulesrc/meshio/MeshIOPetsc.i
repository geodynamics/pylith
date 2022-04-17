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
 * @file modulesrc/meshio/MeshIOPetsc.i
 *
 * @brief Python interface to C++ MeshIOPetsc object.
 */

namespace pylith {
    namespace meshio {
        class MeshIOPetsc: public MeshIO
        { // MeshIOPetsc
          // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Constructor
            MeshIOPetsc(void);

            /// Destructor
            ~MeshIOPetsc(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            // PROTECTED METHODS //////////////////////////////////////////////
protected:

            /// Write mesh
            void _write(void) const;

            /// Read mesh
            void _read(void);

        }; // MeshIOPetsc

    } // meshio
} // pylith

// End of file
