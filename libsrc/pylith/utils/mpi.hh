// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/utils/mpi.hh
 *
 * @brief Simple functions for running in parallel.
 */

#if !defined(pylith_utils_mpi_hh)
#define pylith_utils_mpi_hh

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

#endif // pylith_utils_mpi_hh

// End of file
