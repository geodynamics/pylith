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
 * @file modulesrc/meshio/OutputPhysics.i
 *
 * @brief Python interface to C++ OutputPhysics object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::OutputPhysics :
            public pylith::problems::ObserverPhysics,
            public pylith::meshio::OutputObserver {
            // PUBLIC METHODS ///////////////////////////////////////////////
public:

            /// Constructor
            OutputPhysics(void);

            /// Destructor
            virtual ~OutputPhysics(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set names of information fields requested for output.
             *
             * @param[in] names Array of field names.
             * @param[in] numNames Length of array.
             */
            %apply(const char* const* string_list, const int list_len) {
                (const char* names[],
                 const int numNames)
            };
            void setInfoFields(const char* names[],
                               const int numNames);

            %clear(const char* names[], const int numNames);

            /** Get names of information fields requested for output.
             *
             * @returns Array of field names.
             */
            const pylith::string_vector& getInfoFields(void) const;

            /** Set names of data fields requested for output.
             *
             * @param[in] names Array of field names.
             * @param[in] numNames Length of array.
             */
            %apply(const char* const* string_list, const int list_len) {
                (const char* names[],
                 const int numNames)
            };
            void setDataFields(const char* names[],
                               const int numNames);

            %clear(const char* names[], const int numNames);

            /** Get names of data fields requested for output.
             *
             * @returns Array of field names.
             */
            const pylith::string_vector& getDataFields(void) const;

            /** Set time scale.
             *
             * @param[in] value Time scale for dimensionalizing time.
             */
            void setTimeScale(const PylithReal value);

            /** Verify configuration.
             *
             * @param[in] solution Solution field.
             */
            void verifyConfiguration(const pylith::topology::Field& solution) const;

            /** Receive update (subject of observer).
             *
             * @param[in] t Current time.
             * @param[in] tindex Current time step.
             * @param[in] solution Solution at time t.
             * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after
             * initialization).
             */
            void update(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution,
                        const bool infoOnly);

        }; // OutputPhysics

    } // meshio
} // pylith

// End of file
