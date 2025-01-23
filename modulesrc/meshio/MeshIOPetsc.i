// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

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

            /** Set flag for marking Gmsh vertices.
             *
             * @param value True if marking Gmsh vertices.
             */
            void setGmshMarkVertices(const bool value);

            /** Returns true if marking Gmsh vertices, otherwise false.
             *
             * @returns Mesh format.
             */
            bool getGmshMarkVertices(void) const;

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
