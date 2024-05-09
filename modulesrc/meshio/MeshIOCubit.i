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
 * @file modulesrc/meshio/MeshIOCubit.i
 *
 * @brief Python interface to C++ MeshIOCubit object.
 */

namespace pylith {
    namespace meshio {
        class MeshIOCubit: public MeshIO
        { // MeshIOCubit
          // PUBLIC METHODS /////////////////////////////////////////////////
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

            // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

            /// Write mesh
            void _write(void) const;

            /// Read mesh
            void _read(void);

        }; // MeshIOCubit

    } // meshio
} // pylith

// End of file
