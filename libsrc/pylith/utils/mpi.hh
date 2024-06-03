// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include <mpi.h>

namespace pylith {
    namespace utils {
        class MPI {
public:

            /** Is process root (0) process?
             *
             * @returns True on process 0, otherwise false.
             */
            static
            inline
            bool isRoot(void) {
                int rank = 0;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                return rank == 0;
            } // isRoot

        };
    } // utils
} // pylith

// End of file
