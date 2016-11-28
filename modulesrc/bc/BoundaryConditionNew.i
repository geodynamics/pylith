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

        class BoundaryConditionNew
        { // class BoundaryConditionNew

        // PUBLIC METHODS /////////////////////////////////////////////////
public:

        /// Default constructor.
        BoundaryConditionNew(void);

        /// Destructor.
        virtual
        ~BoundaryConditionNew(void);

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

        }; // class BoundaryConditionNew

    } // bc
} // pylith


// End of file
