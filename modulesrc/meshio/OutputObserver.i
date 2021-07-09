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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
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

            /** Set basis order for output.
             *
             * @param[in] value Basis order for output.
             */
            void setOutputBasisOrder(const int value);

            /** Set time scale.
             *
             * @param[in] value Time scale for dimensionalizing time.
             */
            void setTimeScale(const PylithReal value);

        }; // Observer

    } // meshio
} // pylith

// End of file
