// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


namespace pylith {
    namespace petsc {
        class Application {
            // PUBLIC MEMBERS
            // //////////////////////////////////////////////////////////////////////////////////////////////////
          public:
            /** Initialize
             *
             * @param[in] argc Number of command line arguments.
             * @param[in] argv Command line argument.
             */
            static void initialize(int argc,
                                   char** argv);

            /// Finalize.
            static void finalize(void);

            /** Register citations.
             *
             * @param[in] citation BiBTeX citation string.
             */
            static void registerCitation(const char* citation);

            // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
          private:
            Application(void) = delete;
            Application(const Application&) = delete;
            const Application& operator=(const Application&) = delete;

        }; // class PetscDefaults

    } // namespace petsc
} // namespace pylith


// End of file
