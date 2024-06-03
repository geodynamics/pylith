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
 * @file modulesrc/meshio/OutputSoln.i
 *
 * @brief Python interface to C++ OutputSoln object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::OutputSoln:
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

            /** Receive update (subject of observer).
             *
             * @param[in] t Current time.
             * @param[in] tindex Current time step.
             * @param[in] solution Solution at time t.
             * @param[in] notification Type of notification.
             */
            void update(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution,
                        const pylith::problems::Observer::NotificationType notification);

        }; // OutputSoln

    } // meshio
} // pylith

// End of file
