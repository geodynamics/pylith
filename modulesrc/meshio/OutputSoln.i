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
 * @file modulesrc/meshio/OutputSoln.i
 *
 * @brief Python interface to C++ OutputSoln object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::OutputSoln :
            public pylith::problems::ObserverSoln,
            public pylith::meshio::OutputObserver {
            // PUBLIC METHODS ///////////////////////////////////////////////
public:

            /// Constructor
            OutputSoln(void);

            /// Destructor
            virtual ~OutputSoln(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set names of solution subfields requested for output.
             *
             * @param[in] names Array of subfield names.
             * @param[in] numNames Length of array.
             */
            %apply(const char* const* string_list, const int list_len) {
                (const char* names[],
                 const int numNames)
            };
            void setOutputSubfields(const char* names[],
                                    const int numNames);

            %clear(const char* names[], const int numNames);

            /** Get names of solution subfields requested for output.
             *
             * @returns Array of subfield names.
             */
            const pylith::string_vector& getOutputSubfields(void) const;

            /** Set time scale.
             *
             * @param[in] value Time scale for dimensionalizing time.
             */
            void setTimeScale(const PylithReal value);

            /** Verify observer is compatible with solution.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void verifyConfiguration(const pylith::topology::Field& solution) const;

            /** Receive update from subject.
             *
             * @param[in] t Current time.
             * @param[in] tindex Current time step.
             * @param[in] solution Solution at time t.
             * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after
             * initialization).
             */
            void update(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution);

        }; // OutputSoln

    } // meshio
} // pylith

// End of file
