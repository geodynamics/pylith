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
 * @file modulesrc/utils/PetscVersion.i
 *
 * @brief C++ object for PETSc version information.
 */

namespace pylith {
    namespace utils {
        class PetscVersion
        { // PetscVersion
          // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

            /// Default constructor.
            PetscVersion(void);

            /// Default destrictor.
            ~PetscVersion(void);

            /** Is source from a release?
             *
             * @returns True if source code comes from a release?
             */
            static
            bool isRelease(void);

            /** Get version number.
             *
             * @returns Version number.
             */
            static
            const char* version(void);

            /** Get GIT revision.
             *
             * @returns GIT revision.
             */
            static
            const char* gitRevision(void);

            /** Get date of GIT revision.
             *
             * @returns Date of GIT revision.
             */
            static
            const char* gitDate(void);

            /** Get GIT branch.
             *
             * @returns GIT branch.
             */
            static
            const char* gitBranch(void);

            /** Get PETSC_DIR.
             *
             * @returns PETSC_DIR.
             */
            static
            const char* petscDir(void);

            /** Get PETSC_ARCH.
             *
             * @returns PETSC_ARCH.
             */
            static
            const char* petscArch(void);

        }; // PetscVersion

    } // utils
} // pylith

// End of file
