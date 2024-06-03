// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file modulesrc/topology/Jacobian.i
 *
 * @brief Python interface to C++ Jacobian object.
 */

namespace pylith {
    namespace topology {
        class Jacobian
        { // Jacobian
          // PUBLIC MEMBERS /////////////////////////////////////////////////
public:

            /** Default constructor.
             *
             * @param field Field associated with mesh and solution of the problem.
             * @param matrixType Type of PETSc sparse matrix.
             * @param blockOkay True if okay to use block size equal to fiberDim
             * (all or none of the DOF at each point are constrained).
             */
            Jacobian(const Field& field,
                     const char* matrixType = "aij",
                     const bool blockOkay = false);

            /// Destructor.
            ~Jacobian(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Get PETSc matrix.
             *
             * @returns PETSc sparse matrix.
             */
            const PetscMat* matrix(void) const;

            /** Get PETSc matrix.
             *
             * @returns PETSc sparse matrix.
             */
            PetscMat* matrix(void);

            /** Get matrix type.
             *
             * @returns Matrix type.
             */
            const char* matrixType(void) const;

            /** Assemble matrix.
             *
             * @param mode Assembly mode.
             */
            void assemble(const char* mode);

            /// Set entries in matrix to zero (retain structure).
            void zero(void);

            /// View matrix to stdout.
            void view(void) const;

            /** Write matrix to binary file.
             *
             * @param filename Name of file.
             * @param comm MPI communicator.
             */
            void write(const char* filename,
                       const MPI_Comm comm);

            /// Verify symmetry of matrix. For debugger purposes only.
            void verifySymmetry(void) const;

        }; // Jacobian

    } // topology
} // pylith

// End of file
