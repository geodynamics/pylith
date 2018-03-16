// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/bc/BoundaryCondition.i
 *
 * @brief Python interface to C++ BoundaryCondition object.
 */

namespace pylith {
    namespace bc {

        class BoundaryCondition { // class BoundaryCondition

            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Default constructor.
            BoundaryCondition(void);

            /// Destructor.
            virtual ~BoundaryCondition(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set label of boundary condition surface.
             *
             * @param value Label of surface (from mesh generator).
             */
            void label(const char* value);

            /** Get label of boundary condition surface.
             *
             * @returns Label of surface (from mesh generator).
             */
            const char* label(void) const;

            /** Set name of field in solution to constrain.
             *
             * @param[in] value Name of field in solution to constrain.
             */
            void field(const char* value);

            /** Get name of field in solution to constrain.
             *
             * @returns Name of field in solution to constrain.
             */
            const char* field(void) const;

        }; // class BoundaryCondition

    } // bc
} // pylith


// End of file
