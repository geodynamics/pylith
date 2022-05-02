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
