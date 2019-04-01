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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/meshio/OutputObserver.i
 *
 * @brief Python interface to C++ OutputObserver object.
 */

namespace pylith {
    namespace meshio {
        class OutputObserver : public pylith::utils::PyreComponent {
            // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

            /// Constructor
            OutputObserver(void);

            /// Destructor
            virtual ~OutputObserver(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set output trigger for how often to write output.
             *
             * @param[in] trigger Output trigger.
             */
            void setTrigger(pylith::meshio::OutputTrigger* const trigger);

            /** Get trigger for how often to write otuput.
             *
             * @returns Output trigger.
             */
            const pylith::meshio::OutputTrigger* getTrigger(void) const;

            /** Set writer to write data to file.
             *
             * @param[in] datawriter Writer for data.
             */
            void setWriter(pylith::meshio::DataWriter* const writer);

            /** Set filter for vertex data.
             *
             * @param[in] filter Filter to apply to vertex data before writing.
             */
            void setFieldFilter(pylith::meshio::FieldFilter* const filter);

        }; // Observer

    } // meshio
} // pylith

// End of file
