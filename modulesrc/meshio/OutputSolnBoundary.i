// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/meshio/OutputSolnBoundary.i
 *
 * @brief Python interface to C++ OutputSolnBoundary object.
 */

namespace pylith {
    namespace meshio {

        class pylith::meshio::OutputSolnBoundary : public pylith::meshio::OutputSoln {

            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /** Constructor
             *
             * @param[in] problem Problem to observe.
             */
            OutputSolnBoundary(pylith::problems::Problem* const problem);

            /// Destructor
            virtual ~OutputSolnBoundary(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set label identifier for subdomain.
             *
             * @param[in] value Label of subdomain.
             */
            void label(const char* value);

            /** Verify configuration.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void verifyConfiguration(const pylith::topology::Field& solution) const;

            // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

            /** Write solution at time step.
             *
             * @param[in] t Current time.
             * @param[in] tindex Current time step.
             * @param[in] solution Solution at time t.
             */
            virtual
            void _writeDataStep(const PylithReal t,
                                const PylithInt tindex,
                                const pylith::topology::Field& solution);

        }; // OutputSolnBoundary

    } // meshio
} // pylith


// End of file
