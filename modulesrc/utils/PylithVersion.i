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
 * @file modulesrc/utils/PylithVersion.i
 *
 * @brief C++ object for PyLith version information.
 */

namespace pylith {
    namespace utils {
        class PylithVersion
        { // PylithVersion
          // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

            /// Default constructor.
            PylithVersion(void);

            /// Default destrictor.
            ~PylithVersion(void);

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

            /** Get DOI.
             *
             * @returns DOI.
             */
            static
            const char* doi(void);

            /** Get GIT revision.
             *
             * @returns GIT revision.
             */
            static
            const char* gitRevision(void);

            /** Get GIT hash.
             *
             * @returns GIT hash.
             */
            static
            const char* gitHash(void);

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

        }; // PylithVersion

    } // utils
} // pylith

// End of file
