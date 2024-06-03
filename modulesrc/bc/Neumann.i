// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
/** @file modulesrc/bc/NeumannNew.i
 *
 * @brief Python interface to C++ Neumann object.
 */

namespace pylith {
    namespace bc {
        class Neumann:
        public pylith::feassemble::IntegratorBoundary {
            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Default constructor.
            Neumann(void);

            /// Destructor.
            virtual
            ~Neumann(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Name of scale associated with Neumann boundary
             * condition (e.g., 'pressure' for elasticity).
             *
             * A Neumann boundary condition constrains the gradient in
             * a solution subfield. In some cases the constraint is
             * actually on a scaled version of the gradient as is the
             * case of a Neumann boundary condition for elasticity
             * that constrains boundary tractions.
             *
             * @param value Name of scale for nondimensionalizing Neumann boundary condition.
             */
            void scaleName(const char* value);

        }; // class Neumann

    } // bc
} // pylith

// End of file
